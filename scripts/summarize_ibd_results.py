#!/usr/bin/env python3
"""Summarize parent-child IBD results from GERMLINE and Beagle across chromosomes.

Produces (in --outdir):
    parent_child_check.tsv   - per-chromosome metrics for the specified pair
"""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


# hg19 chromosome lengths (bp) for autosomes
HG19_CHR_LENGTHS = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
    6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
    11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
    16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51244237,
}


def norm_pair(a, b):
    return (a, b) if a <= b else (b, a)


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def intersect_len(a, b):
    a = merge_intervals(a)
    b = merge_intervals(b)
    i = j = total = 0
    while i < len(a) and j < len(b):
        s1, e1 = a[i]
        s2, e2 = b[j]
        s, e = max(s1, s2), min(e1, e2)
        if s < e:
            total += e - s
        if e1 <= e2:
            i += 1
        else:
            j += 1
    return total


def total_bp(intervals):
    return sum(e - s for s, e in merge_intervals(intervals))


def _extract_iid(hap_id):
    """Strip haplotype suffixes used by GERMLINE/germline2."""
    if hap_id.endswith(".0") or hap_id.endswith(".1"):
        return hap_id[:-2]
    if "_" in hap_id:
        parts = hap_id.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]
    return hap_id


def parse_germline(path):
    pairs = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")

            # germline2 format: ID1 ID2 P0 P1 cM #words #gaps
            if len(cols) >= 4 and cols[2].isdigit() and cols[3].isdigit():
                iid1 = _extract_iid(cols[0])
                iid2 = _extract_iid(cols[1])
                start = int(cols[2])
                end = int(cols[3])
                pairs[norm_pair(iid1, iid2)].append((start, end))
                continue

            # legacy GERMLINE format fallback
            p1 = cols[0].split()
            p2 = cols[1].split() if len(cols) > 1 else []
            if len(p1) < 2 or len(p2) < 2 or len(cols) < 4:
                continue
            iid1 = _extract_iid(p1[1])
            iid2 = _extract_iid(p2[1])
            se = cols[3].split()
            if len(se) < 2:
                continue
            pairs[norm_pair(iid1, iid2)].append((int(se[0]), int(se[1])))
    return pairs


def parse_beagle(path):
    pairs = defaultdict(list)
    with gzip.open(path, "rb") as f:
        for line in f:
            line = line.decode().strip()
            if not line:
                continue
            sid1, _h1, sid2, _h2, _chrom, start, end, _lod = line.split("\t")
            pairs[norm_pair(sid1, sid2)].append((int(start), int(end)))
    return pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chromosomes", nargs="+", required=True, help="Chromosome numbers (e.g. 21 22)")
    p.add_argument("--germline-matches", nargs="+", required=True, help="GERMLINE .match files, one per chromosome")
    p.add_argument("--beagle-ibds", nargs="+", required=True, help="Beagle .ibd.gz files, one per chromosome")
    p.add_argument("--parent", default="NA20317")
    p.add_argument("--child", default="NA20318")
    p.add_argument("--outdir", default="results/summary")
    args = p.parse_args()

    if len(args.chromosomes) != len(args.germline_matches) or len(args.chromosomes) != len(args.beagle_ibds):
        raise ValueError("--chromosomes, --germline-matches, and --beagle-ibds must have the same number of entries")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pc_pair = norm_pair(args.parent, args.child)

    pc_rows = []

    for chrom, germ_path, beagle_path in zip(args.chromosomes, args.germline_matches, args.beagle_ibds):
        print(f"Chr{chrom}: parsing GERMLINE ({germ_path})...")
        germ = parse_germline(germ_path)
        print(f"Chr{chrom}: parsing Beagle ({beagle_path})...")
        beagle = parse_beagle(beagle_path)

        chrom_int = int(chrom)
        if chrom_int not in HG19_CHR_LENGTHS:
            raise ValueError(f"Chromosome length not available for chr{chrom} in hg19 table")
        chromosome_bp = HG19_CHR_LENGTHS[chrom_int]

        g_pc = germ.get(pc_pair, [])
        b_pc = beagle.get(pc_pair, [])
        germline_chr_covered_bp = total_bp(g_pc)
        beagle_chr_covered_bp = total_bp(b_pc)
        germline_chr_covered_pct = round(100 * germline_chr_covered_bp / chromosome_bp, 6)
        beagle_chr_covered_pct = round(100 * beagle_chr_covered_bp / chromosome_bp, 6)
        overlap_bp = intersect_len(g_pc, b_pc)
        overlap_chr_pct = round(100 * overlap_bp / chromosome_bp, 6)
        union_bp = germline_chr_covered_bp + beagle_chr_covered_bp - overlap_bp
        jaccard = round(overlap_bp / union_bp, 6) if union_bp > 0 else 0.0

        pc_rows.append({
            "chr": chrom,
            "parent": args.parent,
            "child": args.child,
            "chromosome_bp": chromosome_bp,
            "germline_segments": len(g_pc),
            "beagle_segments": len(b_pc),
            "germline_bp": germline_chr_covered_bp,
            "beagle_bp": beagle_chr_covered_bp,
            "germline_chr_covered_pct": germline_chr_covered_pct,
            "beagle_chr_covered_pct": beagle_chr_covered_pct,
            "overlap_bp": overlap_bp,
            "overlap_chr_pct": overlap_chr_pct,
            "jaccard": jaccard,
        })

    pc_fields = ["chr", "parent", "child", "chromosome_bp", "germline_segments", "beagle_segments",
                 "germline_bp", "beagle_bp", "germline_chr_covered_pct", "beagle_chr_covered_pct",
                 "overlap_bp", "overlap_chr_pct", "jaccard"]

    autosome_rows = []
    for row in pc_rows:
        try:
            chrom_int = int(row["chr"])
        except (TypeError, ValueError):
            continue
        if 1 <= chrom_int <= 22:
            autosome_rows.append(row)

    if autosome_rows:
        avg_row = {
            "chr": "autosome_avg",
            "parent": args.parent,
            "child": args.child,
            "chromosome_bp": "",
            "germline_segments": "",
            "beagle_segments": "",
            "germline_bp": "",
            "beagle_bp": "",
            "germline_chr_covered_pct": round(
                sum(r["germline_chr_covered_pct"] for r in autosome_rows) / len(autosome_rows), 6
            ),
            "beagle_chr_covered_pct": round(
                sum(r["beagle_chr_covered_pct"] for r in autosome_rows) / len(autosome_rows), 6
            ),
            "overlap_bp": "",
            "overlap_chr_pct": round(
                sum(r["overlap_chr_pct"] for r in autosome_rows) / len(autosome_rows), 6
            ),
            "jaccard": round(
                sum(r["jaccard"] for r in autosome_rows) / len(autosome_rows), 6
            ),
        }
        pc_rows.append(avg_row)

    with open(outdir / "parent_child_check.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pc_fields, delimiter="\t")
        w.writeheader()
        w.writerows(pc_rows)

    print(f"Wrote: {outdir}/parent_child_check.tsv")


if __name__ == "__main__":
    main()
