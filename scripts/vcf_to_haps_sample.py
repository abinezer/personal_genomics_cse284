#!/usr/bin/env python3
"""Convert a phased VCF to SHAPEIT/IMPUTE .haps and .sample files.

This prepares germline2-compatible inputs:
  - <out>.haps
  - <out>.sample

Usage:
    python3 vcf_to_haps_sample.py input.vcf.gz --out prefix
"""

import argparse
import gzip
import sys
from pathlib import Path


def open_vcf(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(gt: str):
    if "|" in gt:
        a0s, a1s = gt.split("|", 1)
    elif "/" in gt:
        a0s, a1s = gt.split("/", 1)
    else:
        a0s = a1s = gt

    if a0s in {"0", "1"} and a1s in {"0", "1"}:
        return a0s, a1s
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="Input phased VCF (.vcf or .vcf.gz)")
    parser.add_argument("--out", required=True, help="Output prefix")
    args = parser.parse_args()

    out = Path(args.out)
    haps_path = Path(str(out) + ".haps")
    sample_path = Path(str(out) + ".sample")

    samples = []
    n_written = 0
    skipped = 0

    print("Reading VCF and writing .haps...", file=sys.stderr, flush=True)
    with open_vcf(args.vcf) as fin, open(haps_path, "w") as haps:
        for line in fin:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                continue

            cols = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]

            if "," in alt or len(ref) != 1 or len(alt) != 1:
                skipped += 1
                continue

            fmt = cols[8].split(":")
            if "GT" not in fmt:
                raise ValueError("VCF FORMAT does not contain GT")
            gt_idx = fmt.index("GT")

            hap_bits = []
            valid = True
            for sample_field in cols[9:]:
                gt = sample_field.split(":")[gt_idx]
                parsed = parse_gt(gt)
                if parsed is None:
                    valid = False
                    break
                hap_bits.extend(parsed)

            if not valid:
                skipped += 1
                continue

            marker_id = vid if vid != "." else f"{chrom}:{pos}:{ref}:{alt}"
            haps.write(f"{chrom} {marker_id} {pos} {ref} {alt} {' '.join(hap_bits)}\n")
            n_written += 1

    print("Writing .sample...", file=sys.stderr, flush=True)
    with open(sample_path, "w") as sample:
        sample.write("ID_1 ID_2 missing\n")
        sample.write("0 0 0\n")
        for sid in samples:
            sample.write(f"{sid} {sid} 0\n")

    print(f"Done. Wrote {haps_path} and {sample_path}", file=sys.stderr, flush=True)
    print(f"  Samples: {len(samples)}", file=sys.stderr, flush=True)
    print(f"  Variants written: {n_written}", file=sys.stderr, flush=True)
    if skipped:
        print(f"  Variants skipped (non-biallelic/non-SNP/missing GT): {skipped}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
