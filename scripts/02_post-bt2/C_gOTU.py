#!/usr/bin/env python3
import os
import sys
from biom import load_table
from biom.util import biom_open

# Default Wol taxonomy map
DEFAULT_MAP = "/projects/wol/qiyun/wol2/taxonomy/taxid.map"

def load_taxid_to_gid(map_fp):
    """
    Input lines: GID<TAB>TaxID
    Returns: dict {TaxID: first-seen GID}
    """
    tax2gid = {}
    dup = {}
    with open(map_fp, "r") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            gid, taxid = parts[0].strip(), parts[1].strip()
            if taxid not in tax2gid:
                tax2gid[taxid] = gid
            else:
                dup[taxid] = dup.get(taxid, 1) + 1
    if dup:
        sys.stderr.write(f"[WARN] {len(dup)} TaxIDs map to multiple GIDs; using first seen per TaxID.\n")
    return tax2gid

def remap_biom(in_fp, out_fp, tax2gid):
    tbl = load_table(in_fp)

    # Build explicit id_map: old_id -> new_id (default to old if missing)
    id_map = {}
    missing = 0
    for old in tbl.ids(axis='observation'):
        old_s = old.decode() if isinstance(old, bytes) else str(old)
        new = tax2gid.get(old_s)
        if new is None:
            id_map[old_s] = old_s
            missing += 1
        else:
            id_map[old_s] = new

    # Update IDs in-place
    tbl.update_ids(id_map, axis='observation', inplace=True)

    # Write out
    with biom_open(out_fp, 'w') as h5f:
        tbl.to_hdf5(h5f, generated_by="C_gOTU.py")

    if missing:
        sys.stderr.write(f"[INFO] {missing} features in {os.path.basename(in_fp)} not found in map; left unchanged.\n")

def main():
    # Same TBL_DIRS as before
    TBL_BASE = "/ddn_scratch/yac027/Bowtie2Benchmark/tables/Test1"
    TBL_DIRS = [
        "01_k16_no-split",
        "02_a_no-split",
        "03_igor_no-split",
        "04_k16_split4",
        "05_a_split4",
        "06_igor_split4",
        "07_k16_split10",
        "08_a_split10",
        "09_igor_split10",
        "04a_k16_split4a",
        "05a_a_split4a",
        "06a_igor_split4a",
    ]

    tax2gid = load_taxid_to_gid(DEFAULT_MAP)

    for d in TBL_DIRS:
        tbl_dir = os.path.join(TBL_BASE, d)

        for level in ["phylum", "genus", "species"]:
            in_fp = os.path.join(tbl_dir, f"{level}_collapsed.biom")
            out_fp = os.path.join(tbl_dir, f"{level}_collapsed_gOTU.biom")

            if not os.path.isfile(in_fp):
                sys.stderr.write(f"[WARN] {in_fp} not found, skipping.\n")
                continue

            sys.stderr.write(f"[INFO] Remapping {in_fp} -> {out_fp}\n")
            remap_biom(in_fp, out_fp, tax2gid)

    sys.stderr.write("[DONE] Finished gOTU remapping for all collapsed tables.\n")

if __name__ == "__main__":
    main()

