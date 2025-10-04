#!/usr/bin/env python3
import os
import sys
import re
from biom import load_table

# Base path and dirs to process
TBL_BASE = "/ddn_scratch/yac027/03_Bowtie2Benchmark/tables/Test1"
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
    "06a_igor_split4a"
]

# Pattern to strip index suffix (_index1, _index2, etc.)
def strip_index(sample_id: str) -> str:
    return re.sub(r"_index\d+$", "", sample_id)

def collapse_samples(in_fp, out_fp, mapping):
    tbl = load_table(in_fp)

    collapsed = tbl.collapse(
        lambda id_, md: mapping[id_],
        axis="sample",
        norm=False
    )

    with open(out_fp, "w") as f:
        collapsed.to_json("collapsed", f)

    return tbl, collapsed   # return both for logging

def main():
    biom_files = ["phylum.biom", "genus.biom", "species.biom"]

    for d in TBL_DIRS:
        dir_path = os.path.join(TBL_BASE, d)
        if not os.path.isdir(dir_path):
            sys.stderr.write(f"[WARN] Missing directory {dir_path}, skipping.\n")
            continue

        for biom_file in biom_files:
            in_fp = os.path.join(dir_path, biom_file)
            if not os.path.exists(in_fp):
                sys.stderr.write(f"[WARN] Missing file {in_fp}, skipping.\n")
                continue

            base, ext = os.path.splitext(in_fp)
            out_fp = f"{base}_collapsed{ext or '.biom'}"

            # Build mapping: strip _indexN suffix → collapsed sample ID
            tbl = load_table(in_fp)
            mapping = {sid: strip_index(sid) for sid in tbl.ids(axis="sample")}

            # Collapse and log shapes
            orig_tbl, collapsed_tbl = collapse_samples(in_fp, out_fp, mapping)

            sys.stderr.write(
                f"[INFO] Collapsing {in_fp} "
                f"(input {orig_tbl.shape[0]} features x {orig_tbl.shape[1]} samples) "
                f"-> (output {collapsed_tbl.shape[0]} features x {collapsed_tbl.shape[1]} samples)\n"
            )

    sys.stderr.write("[DONE] Finished collapsing all tables.\n")

if __name__ == "__main__":
    main()

