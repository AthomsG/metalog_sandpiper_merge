#!/usr/bin/env python3
"""
Fetch metalog metadata + the latest sandpiper profiles.

Metalog files (small, always re-downloaded) -> fetched_data/:
  - human_extended_wide_latest.tsv.gz   (human environment metadata, wide)
  - ocean_extended_wide_latest.tsv.gz   (ocean environment metadata, wide)
  - sequencing_db_mapping_latest.tsv.gz (metalog sample_alias <-> run accession map)

Sandpiper profiles (large, ~3.3 GB) -> SANDPIPER_DIR:
  - Resolves the LATEST version from the Zenodo concept record and downloads the
    GTDB profile file (sandpiper<version>.gtdb.csv.gz).
  - The version is part of the filename, so a download only happens when a *newer*
    version is published (the current versioned file is otherwise found locally and
    skipped, after a size check).
  - Verifies the MD5 of freshly downloaded files.
  - Writes the resolved path to fetched_data/sandpiper_latest.txt, which
    02_merge_and_count.py reads to locate the sandpiper file.
"""

import hashlib
import json
import os
import urllib.request

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FETCHED_DIR = os.path.join(HERE, "fetched_data")

# --- metalog ---
BASE = "https://metalog.embl.de/static/download"
METALOG_DOWNLOADS = {
    f"{BASE}/metadata/human_extended_wide_latest.tsv.gz": "human_extended_wide_latest.tsv.gz",
    f"{BASE}/metadata/ocean_extended_wide_latest.tsv.gz": "ocean_extended_wide_latest.tsv.gz",
    f"{BASE}/sequencing_db_mapping_latest.tsv.gz": "sequencing_db_mapping_latest.tsv.gz",
}

# --- sandpiper ---
# Zenodo "latest version" endpoint for the sandpiper concept (concept record 10547493).
# Querying any version's /versions/latest resolves to the newest published version.
SANDPIPER_LATEST_API = "https://zenodo.org/api/records/20437114/versions/latest"
# Which profile file to grab from the record. GTDB taxonomy is the current default;
# change to e.g. "globdb.csv.gz" to use a different taxonomy.
SANDPIPER_PATTERN = "gtdb.csv.gz"
# Where the (large) sandpiper file is stored locally.
SANDPIPER_DIR = "/lisc/data/work/dome/pollak/gaehtgens/data/sandpiper"
SANDPIPER_POINTER = os.path.join(FETCHED_DIR, "sandpiper_latest.txt")
VERIFY_MD5 = True


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------
def _human_size(n) -> str:
    f = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if f < 1024 or unit == "GB":
            return f"{f:.1f} {unit}"
        f /= 1024


def _download(url: str, dest: str):
    # Report progress on its own line every ~10% (no carriage-return spam in logs).
    state = {"next_pct": 0}

    def hook(blocks, bs, total):
        if total > 0:
            done = min(blocks * bs, total)
            pct = done * 100 // total
            if pct >= state["next_pct"]:
                print(f"    {pct:3d}%  {_human_size(done)} / {_human_size(total)}", flush=True)
                state["next_pct"] = pct - (pct % 10) + 10

    urllib.request.urlretrieve(url, dest, reporthook=hook)


def _md5(path: str, chunk: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# ------------------------------------------------------------
# METALOG
# ------------------------------------------------------------
def fetch_metalog():
    for url, fname in METALOG_DOWNLOADS.items():
        dest = os.path.join(FETCHED_DIR, fname)
        print(f"Downloading {url}\n        -> {dest}", flush=True)
        urllib.request.urlretrieve(url, dest)  # always re-download to stay current
        print(f"  done: {_human_size(os.path.getsize(dest))}", flush=True)


# ------------------------------------------------------------
# SANDPIPER
# ------------------------------------------------------------
def fetch_sandpiper():
    print("Resolving latest sandpiper version from Zenodo ...", flush=True)
    with urllib.request.urlopen(SANDPIPER_LATEST_API) as r:
        meta = json.load(r)

    version = meta.get("metadata", {}).get("version", "?")
    record_id = meta["id"]

    matches = [f for f in meta["files"] if f["key"].endswith(SANDPIPER_PATTERN)]
    if not matches:
        raise RuntimeError(
            f"No file matching '*{SANDPIPER_PATTERN}' in latest sandpiper record {record_id}"
        )
    f = matches[0]
    key = f["key"]
    size = f["size"]
    md5_expected = f["checksum"].split(":", 1)[-1]  # 'md5:...' -> '...'
    url = f"https://zenodo.org/records/{record_id}/files/{key}?download=1"

    os.makedirs(SANDPIPER_DIR, exist_ok=True)
    dest = os.path.join(SANDPIPER_DIR, key)

    print(f"Latest sandpiper version: {version}  (record {record_id})", flush=True)
    print(f"Target file: {key}  ({_human_size(size)})", flush=True)

    if os.path.exists(dest) and os.path.getsize(dest) == size:
        print(f"  already present locally (size matches), skipping download:\n    {dest}", flush=True)
    else:
        print(f"  downloading -> {dest}", flush=True)
        _download(url, dest)
        actual = os.path.getsize(dest)
        if actual != size:
            raise RuntimeError(f"size mismatch for {key}: got {actual}, expected {size}")
        if VERIFY_MD5:
            print("  verifying md5 ...", flush=True)
            actual_md5 = _md5(dest)
            if actual_md5 != md5_expected:
                raise RuntimeError(
                    f"md5 mismatch for {key}: got {actual_md5}, expected {md5_expected}"
                )
            print("  md5 OK", flush=True)

    with open(SANDPIPER_POINTER, "w") as out:
        out.write(dest + "\n")
    print(f"  pointer written -> {SANDPIPER_POINTER}", flush=True)


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main():
    os.makedirs(FETCHED_DIR, exist_ok=True)
    fetch_metalog()
    fetch_sandpiper()
    print("All downloads complete.", flush=True)


if __name__ == "__main__":
    main()
