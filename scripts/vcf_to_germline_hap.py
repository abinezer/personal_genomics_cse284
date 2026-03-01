#!/usr/bin/env python3
"""Convert a phased VCF to GERMLINE haploid PED/MAP format.

Each individual in the VCF becomes TWO rows in the output PED file,
one per haplotype. For GERMLINE's haploid diploid-trick, each allele
is written twice (GERMLINE expects allele pairs even in haploid mode).

Usage:
    python3 vcf_to_germline_hap.py input.vcf.gz --out prefix
"""
import argparse
import gzip
import sys
from pathlib import Path


def open_vcf(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("vcf", help="Phased VCF file (.vcf or .vcf.gz)")
    p.add_argument("--out", required=True, help="Output prefix (produces .ped and .map)")
    args = p.parse_args()

    samples: list[str] = []
    var_info: list[tuple[str, str, str]] = []  # (chrom, id, pos)
    ref_alleles: list[str] = []
    alt_alleles: list[str] = []
    # hap0[variant_idx] = [allele_int per sample]  (0=ref, 1=alt, -1=missing)
    hap0: list[list[int]] = []
    hap1: list[list[int]] = []
    skipped_non_snp_or_multiallelic = 0

    print("Reading VCF...", file=sys.stderr, flush=True)
    with open_vcf(args.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                print(f"  {len(samples)} samples", file=sys.stderr, flush=True)
                continue

            cols = line.rstrip("\n").split("\t")
            chrom = cols[0]
            pos = cols[1]
            rsid = cols[2] if cols[2] != "." else f"chr{chrom}:{pos}"
            ref = cols[3]
            alt = cols[4]

            # Only keep biallelic SNPs. Anything else is unsafe for this
            # converter and should be excluded upstream or here.
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                skipped_non_snp_or_multiallelic += 1
                continue

            var_info.append((chrom, rsid, pos))
            ref_alleles.append(ref)
            alt_alleles.append(alt)

            fmt = cols[8].split(":")
            if "GT" not in fmt:
                raise ValueError("VCF FORMAT field does not contain GT")
            gt_idx = fmt.index("GT")
            h0_row: list[int] = []
            h1_row: list[int] = []
            for gts in cols[9:]:
                gt = gts.split(":")[gt_idx]
                if "|" in gt:
                    a0s, a1s = gt.split("|", 1)
                elif "/" in gt:
                    a0s, a1s = gt.split("/", 1)
                else:
                    a0s = a1s = gt
                a0 = int(a0s) if a0s.isdigit() else -1
                a1 = int(a1s) if a1s.isdigit() else -1
                if a0 not in (-1, 0, 1) or a1 not in (-1, 0, 1):
                    # Multi-allelic genotype slipped through; skip variant.
                    h0_row = []
                    h1_row = []
                    skipped_non_snp_or_multiallelic += 1
                    break
                h0_row.append(a0)
                h1_row.append(a1)
            if not h0_row:
                continue
            hap0.append(h0_row)
            hap1.append(h1_row)

    n_var = len(var_info)
    n_sam = len(samples)
    print(f"  {n_var} variants", file=sys.stderr, flush=True)
    if skipped_non_snp_or_multiallelic:
        print(
            f"  skipped {skipped_non_snp_or_multiallelic} non-biallelic/non-SNP variants",
            file=sys.stderr,
            flush=True,
        )

    out = Path(args.out)

    # --- Write MAP file ---
    print("Writing .map...", file=sys.stderr, flush=True)
    with open(str(out) + ".map", "w") as f:
        for chrom, rsid, pos in var_info:
            f.write(f"{chrom}\t{rsid}\t0\t{pos}\n")

    # --- Write haploid PED file ---
    # For each individual, emit TWO rows (one per haplotype).
    # Each marker column is duplicated so GERMLINE sees homozygous markers
    # for that haplotype (haploid-diploid trick).
    print("Writing haploid .ped...", file=sys.stderr, flush=True)
    # Pre-transpose: build per-sample haplotype allele lists
    # hap0_T[sample_idx] = list of allele ints across variants
    hap0_T: list[list[int]] = [[] for _ in range(n_sam)]
    hap1_T: list[list[int]] = [[] for _ in range(n_sam)]
    for var_idx in range(n_var):
        h0 = hap0[var_idx]
        h1 = hap1[var_idx]
        for s in range(n_sam):
            hap0_T[s].append(h0[s])
            hap1_T[s].append(h1[s])

    with open(str(out) + ".ped", "w") as f:
        for s_idx, sid in enumerate(samples):
            for hap_label, hap_alleles in enumerate([hap0_T[s_idx], hap1_T[s_idx]], start=1):
                row_parts = [sid, f"{sid}_{hap_label}", "0", "0", "0", "-9"]
                for a_int, ref, alt in zip(hap_alleles, ref_alleles, alt_alleles):
                    if a_int == 0:
                        allele = ref
                    elif a_int == 1:
                        allele = alt
                    else:
                        allele = "0"
                    # Duplicate allele for GERMLINE's diploid expectations
                    row_parts.append(allele)
                    row_parts.append(allele)
                f.write("\t".join(row_parts) + "\n")
            if (s_idx + 1) % 10 == 0:
                print(f"  {s_idx + 1}/{n_sam} samples written", file=sys.stderr, flush=True)

    print(f"Done. Output: {out}.ped, {out}.map", file=sys.stderr, flush=True)
    print(f"  Haplotypes: {2 * n_sam}, Variants: {n_var}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
