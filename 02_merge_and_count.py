#!/usr/bin/env python3
"""
Per-environment (human / ocean) shared runs + species occurrences.

Join chain:
    wide.sample_alias  ->  mapping.sample_alias (kind=="run")  ->  mapping.external_id
    (run accession)    ->  sandpiper.sample

For each environment we output, into output/<env>/:
  - shared_runs.txt          : run accessions present in BOTH that env's metalog set
                               and sandpiper (one per line, sorted).
  - species_occurrences.csv  : columns [species, occurrences]; occurrences = number of
                               sandpiper rows (one per species per run) in that env's
                               shared runs where the s__ species appears.

Sandpiper is parsed exactly once for both environments.
"""

import glob
import gzip
import os
import shutil
import subprocess
from collections import defaultdict

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FETCHED_DIR = os.path.join(HERE, "fetched_data")
OUT_DIR = os.path.join(HERE, "output")

MAPPING_FILE = os.path.join(FETCHED_DIR, "sequencing_db_mapping_latest.tsv.gz")

# Sandpiper file is downloaded by 01_fetch_data.py. Its path is recorded in this pointer
# file; fall back to globbing the sandpiper data dir for the latest GTDB profile.
SANDPIPER_POINTER = os.path.join(FETCHED_DIR, "sandpiper_latest.txt")
SANDPIPER_DIR = "/lisc/data/work/dome/pollak/gaehtgens/data/sandpiper"
SANDPIPER_GLOB = "*.gtdb.csv.gz"

# env name -> wide metadata file
ENV_WIDE_FILES = {
    "human": os.path.join(FETCHED_DIR, "human_extended_wide_latest.tsv.gz"),
    "ocean": os.path.join(FETCHED_DIR, "ocean_extended_wide_latest.tsv.gz"),
}

USE_PIGZ_IF_AVAILABLE = True
PROGRESS_EVERY_LINES = 1_000_000

# Only count a species in a sandpiper row when that row's coverage is >= this value.
# Does NOT affect which runs are listed as shared (only the species counts).
MIN_COVERAGE = 10.0


# ------------------------------------------------------------
# IO HELPERS
# ------------------------------------------------------------
def open_maybe_pigz(path: str):
    """Return (proc, fh). Uses pigz for faster decompression when available."""
    if USE_PIGZ_IF_AVAILABLE and shutil.which("pigz"):
        proc = subprocess.Popen(
            ["pigz", "-dc", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024,
        )
        return proc, proc.stdout
    return None, gzip.open(path, "rt")


def _finish_reader(proc, fh):
    if fh is not None:
        fh.close()
    if proc is not None:
        _, stderr = proc.communicate()
        if proc.returncode not in (0, None):
            raise RuntimeError(f"pigz failed ({proc.returncode}): {stderr}")


def extract_species_from_taxonomy(taxonomy: str):
    """All taxonomy tokens beginning with 's__', kept verbatim; unique within a row."""
    species_found = set()
    for token in taxonomy.split(";"):
        token = token.strip()
        if token.startswith("s__"):
            species_found.add(token)
    return species_found


# ------------------------------------------------------------
# STEP 1: per-environment sample_alias sets from wide metadata
# ------------------------------------------------------------
def load_env_aliases(env_wide_files: dict):
    """env -> set(sample_alias). Locates the sample_alias column by header name."""
    env_aliases = {}
    for env, path in env_wide_files.items():
        aliases = set()
        proc, fh = open_maybe_pigz(path)
        try:
            header = next(fh, None)
            if header is None:
                env_aliases[env] = aliases
                continue
            cols = header.rstrip("\n").split("\t")
            alias_idx = cols.index("sample_alias")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) > alias_idx:
                    a = parts[alias_idx].strip()
                    if a:
                        aliases.add(a)
        finally:
            _finish_reader(proc, fh)
        env_aliases[env] = aliases
        print(f"[aliases] {env}: {len(aliases):,} sample_alias values", flush=True)
    return env_aliases


# ------------------------------------------------------------
# STEP 2: alias -> run accession, restricted to env alias sets
# ------------------------------------------------------------
def build_run_to_envs(mapping_file: str, env_aliases: dict):
    """
    Stream the mapping file once. For each kind=="run" row whose sample_alias is in an
    environment's alias set, attach that env to the run accession (external_id).

    Returns:
        run_to_envs : dict run_accession -> set(env)
        env_runs    : dict env -> set(run_accession)
    """
    run_to_envs = defaultdict(set)
    env_runs = {env: set() for env in env_aliases}

    proc, fh = open_maybe_pigz(mapping_file)
    try:
        header = next(fh, None)
        cols = header.rstrip("\n").split("\t")
        alias_idx = cols.index("sample_alias")
        kind_idx = cols.index("kind")
        ext_idx = cols.index("external_id")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(alias_idx, kind_idx, ext_idx):
                continue
            if parts[kind_idx] != "run":
                continue
            alias = parts[alias_idx].strip()
            run = parts[ext_idx].strip()
            if not run:
                continue
            for env, aliases in env_aliases.items():
                if alias in aliases:
                    run_to_envs[run].add(env)
                    env_runs[env].add(run)
    finally:
        _finish_reader(proc, fh)

    for env in env_aliases:
        print(f"[runs] {env}: {len(env_runs[env]):,} run accessions from metalog", flush=True)
    return run_to_envs, env_runs


# ------------------------------------------------------------
# STEP 3: single sandpiper pass for all environments
# ------------------------------------------------------------
def scan_sandpiper(path: str, run_to_envs: dict, envs):
    """
    One pass over sandpiper. Returns:
      shared_runs        : dict env -> set(run accession seen in sandpiper)
      species_occurrences: dict env -> dict(species -> count)
    """
    shared_runs = {env: set() for env in envs}
    species_occurrences = {env: defaultdict(int) for env in envs}

    proc, fh = open_maybe_pigz(path)
    total_lines = 0
    try:
        next(fh, None)  # header
        for line in fh:
            total_lines += 1

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            sample = parts[0].strip()
            hit_envs = run_to_envs.get(sample)
            if hit_envs:
                # Run membership is independent of coverage.
                for env in hit_envs:
                    shared_runs[env].add(sample)

                # Species counts only from rows with coverage >= MIN_COVERAGE.
                try:
                    coverage = float(parts[1])
                except ValueError:
                    continue
                if coverage < MIN_COVERAGE:
                    continue

                taxonomy = parts[2].strip()
                species_set = extract_species_from_taxonomy(taxonomy) if taxonomy else None
                if species_set:
                    for env in hit_envs:
                        counts = species_occurrences[env]
                        for species in species_set:
                            counts[species] += 1

            if total_lines % PROGRESS_EVERY_LINES == 0:
                summary = "; ".join(
                    f"{env}: runs={len(shared_runs[env]):,}, sp={len(species_occurrences[env]):,}"
                    for env in envs
                )
                print(f"[scan] {total_lines:,} rows | {summary}", flush=True)
    finally:
        _finish_reader(proc, fh)

    summary = "; ".join(
        f"{env}: runs={len(shared_runs[env]):,}, sp={len(species_occurrences[env]):,}"
        for env in envs
    )
    print(f"[scan] done: {total_lines:,} rows | {summary}", flush=True)
    return shared_runs, species_occurrences


# ------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------
def write_env_outputs(env, runs, species_counts):
    env_dir = os.path.join(OUT_DIR, env)
    os.makedirs(env_dir, exist_ok=True)

    runs_path = os.path.join(env_dir, "shared_runs.txt")
    with open(runs_path, "w") as out:
        for run in sorted(runs):
            out.write(run + "\n")
    print(f"[{env}] wrote {len(runs):,} shared runs -> {runs_path}", flush=True)

    csv_path = os.path.join(env_dir, "species_occurrences.csv")
    with open(csv_path, "w") as out:
        out.write("species,occurrences\n")
        for species, count in sorted(species_counts.items(), key=lambda kv: (-kv[1], kv[0])):
            if "," in species or '"' in species:
                species = '"' + species.replace('"', '""') + '"'
            out.write(f"{species},{count}\n")
    print(f"[{env}] wrote {len(species_counts):,} species -> {csv_path}", flush=True)


# ------------------------------------------------------------
# SANDPIPER FILE RESOLUTION
# ------------------------------------------------------------
def resolve_sandpiper_file():
    """Path to the sandpiper file: prefer the pointer written by 01_fetch_data.py,
    else fall back to the newest GTDB profile in SANDPIPER_DIR."""
    if os.path.exists(SANDPIPER_POINTER):
        with open(SANDPIPER_POINTER) as fh:
            path = fh.read().strip()
        if path and os.path.exists(path):
            return path

    candidates = sorted(glob.glob(os.path.join(SANDPIPER_DIR, SANDPIPER_GLOB)))
    if candidates:
        return candidates[-1]

    raise RuntimeError(
        f"No sandpiper file found via {SANDPIPER_POINTER} or {SANDPIPER_DIR}/{SANDPIPER_GLOB}. "
        "Run 01_fetch_data.py first."
    )


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    envs = list(ENV_WIDE_FILES)

    sandpiper_file = resolve_sandpiper_file()
    print(f"MIN_COVERAGE={MIN_COVERAGE} (species counted only from rows with coverage >= this)", flush=True)
    print(f"Sandpiper file: {sandpiper_file}", flush=True)
    print("Loading per-environment sample_alias sets ...", flush=True)
    env_aliases = load_env_aliases(ENV_WIDE_FILES)

    print("Mapping aliases -> run accessions ...", flush=True)
    run_to_envs, env_runs = build_run_to_envs(MAPPING_FILE, env_aliases)

    print("Scanning sandpiper (single pass) ...", flush=True)
    shared_runs, species_occurrences = scan_sandpiper(sandpiper_file, run_to_envs, envs)

    for env in envs:
        write_env_outputs(env, shared_runs[env], species_occurrences[env])

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
