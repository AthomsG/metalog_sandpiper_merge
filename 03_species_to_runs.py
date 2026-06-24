#!/usr/bin/env python3
"""
Given a species name, find the runs (restricted to the human/ocean shared runs used
by 02_merge_and_count.py) where that species is present at coverage >= a threshold.

Usage:
    python3 03_species_to_runs.py "Phocaeicola vulgatus" [--min-coverage 10]

Matching is exact: the species argument (with or without a leading "s__") must match
a sandpiper taxonomy token exactly, e.g. "s__Phocaeicola vulgatus".

Writes, for each environment with at least one match:
    output/species_samples/<env>/<species_slug>.csv   (columns: run,coverage)
"""

import argparse
import glob
import gzip
import os
import re
import shutil
import subprocess
from collections import defaultdict

# ------------------------------------------------------------
# CONFIG (mirrors 02_merge_and_count.py)
# ------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FETCHED_DIR = os.path.join(HERE, "fetched_data")
OUT_DIR = os.path.join(HERE, "output")

MAPPING_FILE = os.path.join(FETCHED_DIR, "sequencing_db_mapping_latest.tsv.gz")

SANDPIPER_POINTER = os.path.join(FETCHED_DIR, "sandpiper_latest.txt")
SANDPIPER_DIR = "/lisc/data/work/dome/pollak/gaehtgens/data/sandpiper"
SANDPIPER_GLOB = "*.gtdb.csv.gz"

ENV_WIDE_FILES = {
    "human": os.path.join(FETCHED_DIR, "human_extended_wide_latest.tsv.gz"),
    "ocean": os.path.join(FETCHED_DIR, "ocean_extended_wide_latest.tsv.gz"),
}

USE_PIGZ_IF_AVAILABLE = True
PROGRESS_EVERY_LINES = 1_000_000

DEFAULT_MIN_COVERAGE = 10.0


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


def resolve_sandpiper_file():
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
# STEP 1: per-environment sample_alias sets, then alias -> run accession
# (same logic as 02_merge_and_count.py)
# ------------------------------------------------------------
def load_env_aliases(env_wide_files: dict):
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


def build_run_to_envs(mapping_file: str, env_aliases: dict):
    run_to_envs = defaultdict(set)

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
    finally:
        _finish_reader(proc, fh)

    return run_to_envs


# ------------------------------------------------------------
# STEP 2: scan sandpiper for the target species
# ------------------------------------------------------------
def find_species_runs(path: str, run_to_envs: dict, target_token: str, min_coverage: float, envs):
    """
    One pass over sandpiper. Returns dict env -> list of (run, coverage) where:
      - run is in run_to_envs (i.e. shared with that env's metalog set)
      - the row's taxonomy contains the exact token target_token
      - coverage >= min_coverage
    """
    results = {env: [] for env in envs}

    proc, fh = open_maybe_pigz(path)
    total_lines = 0
    try:
        next(fh, None)  # header
        for line in fh:
            total_lines += 1

            # Fast pre-filter before splitting/parsing.
            if target_token not in line:
                if total_lines % PROGRESS_EVERY_LINES == 0:
                    print(f"[scan] {total_lines:,} rows processed", flush=True)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            sample = parts[0].strip()
            hit_envs = run_to_envs.get(sample)
            if not hit_envs:
                continue

            try:
                coverage = float(parts[1])
            except ValueError:
                continue
            if coverage < min_coverage:
                continue

            tokens = [t.strip() for t in parts[2].split(";")]
            if target_token not in tokens:
                continue

            for env in hit_envs:
                results[env].append((sample, coverage))

            if total_lines % PROGRESS_EVERY_LINES == 0:
                print(f"[scan] {total_lines:,} rows processed", flush=True)
    finally:
        _finish_reader(proc, fh)

    print(f"[scan] done: {total_lines:,} rows processed", flush=True)
    return results


# ------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------
def slugify(species: str) -> str:
    slug = species.strip().lower().replace(" ", "_")
    return re.sub(r"[^a-z0-9_.-]", "", slug)


def write_species_outputs(species: str, results: dict, min_coverage: float):
    slug = slugify(species)
    out_dir = os.path.join(OUT_DIR, "species_samples")

    # Always write one separate file per environment, even if a list ends up empty.
    for env, rows in results.items():
        env_dir = os.path.join(out_dir, env)
        os.makedirs(env_dir, exist_ok=True)
        csv_path = os.path.join(env_dir, f"{slug}.csv")

        rows_sorted = sorted(rows, key=lambda rc: (-rc[1], rc[0]))
        with open(csv_path, "w") as out:
            out.write("run,coverage\n")
            for run, coverage in rows_sorted:
                out.write(f"{run},{coverage}\n")

        if rows_sorted:
            print(f"[{env}] wrote {len(rows_sorted):,} runs -> {csv_path}", flush=True)
        else:
            print(f"[{env}] no runs found for 's__{species}' at coverage >= {min_coverage} "
                  f"(empty list written -> {csv_path})", flush=True)


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("species", help='Species name, e.g. "Phocaeicola vulgatus" (with or without leading "s__")')
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                         help=f"Minimum coverage required (default: {DEFAULT_MIN_COVERAGE})")
    args = parser.parse_args()

    species = args.species[3:] if args.species.startswith("s__") else args.species
    target_token = f"s__{species}"

    sandpiper_file = resolve_sandpiper_file()
    print(f"Species: {target_token}", flush=True)
    print(f"Min coverage: {args.min_coverage}", flush=True)
    print(f"Sandpiper file: {sandpiper_file}", flush=True)

    envs = list(ENV_WIDE_FILES)

    print("Loading per-environment sample_alias sets ...", flush=True)
    env_aliases = load_env_aliases(ENV_WIDE_FILES)

    print("Mapping aliases -> run accessions ...", flush=True)
    run_to_envs = build_run_to_envs(MAPPING_FILE, env_aliases)

    print("Scanning sandpiper for matching species ...", flush=True)
    results = find_species_runs(sandpiper_file, run_to_envs, target_token, args.min_coverage, envs)

    write_species_outputs(species, results, args.min_coverage)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
