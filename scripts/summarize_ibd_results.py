#!/usr/bin/env python3
"""Summarize and compare IBD results from Beagle, GERMLINE (phased), and PLINK.

Produces:
  overall_metrics.tsv   - aggregate statistics
  pair_metrics.tsv      - per-pair statistics for every pair seen
  parent_child_check.tsv - validation for the known parent-child pair
"""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


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
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        s1, e1 = a[i]
        s2, e2 = b[j]
        s = max(s1, s2)
        e = min(e1, e2)
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
    """Strip haplotype suffix: 'NA20317_1' -> 'NA20317'."""
    if "_" in hap_id:
        parts = hap_id.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]
    return hap_id


def parse_germline_phased(path):
    """Parse GERMLINE match file with haplotype-level IDs.

    Collapses haplotype pairs to individual-level pairs by stripping the
    '_1'/'_2' suffix from IIDs.  Returns dict mapping individual pair ->
    list of (start, end) intervals (union over all haplotype combinations).
    """
    pairs = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            # cols[0] = "FID1 IID1", cols[1] = "FID2 IID2"
            p1 = cols[0].split()
            p2 = cols[1].split()
            if len(p1) < 2 or len(p2) < 2:
                # Skip malformed rows rather than crashing.
                continue
            iid1 = _extract_iid(p1[1])
            iid2 = _extract_iid(p2[1])
            start, end = cols[3].split()
            pairs[norm_pair(iid1, iid2)].append((int(start), int(end)))
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


def parse_plink_genome(path):
    out = {}
    with open(path) as f:
        header = f.readline().strip().split()
        if not header:
            return out
        idx = {h: i for i, h in enumerate(header)}
        required = {"IID1", "IID2", "PI_HAT", "Z0", "Z1", "Z2"}
        if not required.issubset(idx):
            missing = sorted(required - set(idx))
            raise ValueError(f"PLINK genome file missing required columns: {missing}")
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if len(parts) < len(header):
                continue
            key = norm_pair(parts[idx["IID1"]], parts[idx["IID2"]])
            out[key] = {
                "pi_hat": float(parts[idx["PI_HAT"]]),
                "z0": float(parts[idx["Z0"]]),
                "z1": float(parts[idx["Z1"]]),
                "z2": float(parts[idx["Z2"]]),
            }
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beagle-ibd", required=True, help="Beagle .ibd.gz file")
    p.add_argument("--germline-match", required=True, help="GERMLINE .match file (phased haploid)")
    p.add_argument("--plink-genome", required=True, help="PLINK .genome file")
    p.add_argument("--chr-length", type=int, default=51_244_237, help="Chromosome length in bp")
    p.add_argument("--parent", default="NA20317", help="Parent sample ID")
    p.add_argument("--child", default="NA20318", help="Child sample ID")
    p.add_argument("--outdir", default="results/summary", help="Output directory")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    beagle = parse_beagle(Path(args.beagle_ibd))
    germ = parse_germline_phased(Path(args.germline_match))
    plink = parse_plink_genome(Path(args.plink_genome))

    all_pairs = sorted(set(beagle) | set(germ) | set(plink))

    pair_rows = []
    beagle_seg_n = germ_seg_n = 0
    beagle_bp_total = germ_bp_total = overlap_bp_total = 0

    for pair in all_pairs:
        b_int = beagle.get(pair, [])
        g_int = germ.get(pair, [])
        b_bp = total_bp(b_int)
        g_bp = total_bp(g_int)
        o_bp = intersect_len(b_int, g_int)
        pl = plink.get(pair, {})

        beagle_seg_n += len(b_int)
        germ_seg_n += len(g_int)
        beagle_bp_total += b_bp
        germ_bp_total += g_bp
        overlap_bp_total += o_bp

        pair_rows.append({
            "id1": pair[0],
            "id2": pair[1],
            "beagle_segments": len(b_int),
            "germline_segments": len(g_int),
            "beagle_bp": b_bp,
            "germline_bp": g_bp,
            "overlap_bp": o_bp,
            "beagle_chr_prop": round(b_bp / args.chr_length, 6),
            "germline_chr_prop": round(g_bp / args.chr_length, 6),
            "plink_pi_hat": pl.get("pi_hat", ""),
            "plink_z0": pl.get("z0", ""),
            "plink_z1": pl.get("z1", ""),
            "plink_z2": pl.get("z2", ""),
        })

    pair_fields = [
        "id1",
        "id2",
        "beagle_segments",
        "germline_segments",
        "beagle_bp",
        "germline_bp",
        "overlap_bp",
        "beagle_chr_prop",
        "germline_chr_prop",
        "plink_pi_hat",
        "plink_z0",
        "plink_z1",
        "plink_z2",
    ]
    with open(outdir / "pair_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pair_fields, delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    union_bp = beagle_bp_total + germ_bp_total - overlap_bp_total
    overall = {
        "beagle_segment_count": beagle_seg_n,
        "germline_segment_count": germ_seg_n,
        "beagle_total_bp": beagle_bp_total,
        "germline_total_bp": germ_bp_total,
        "pairwise_overlap_bp": overlap_bp_total,
        "jaccard_bp": round(overlap_bp_total / union_bp, 6) if union_bp > 0 else 0.0,
    }
    with open(outdir / "overall_metrics.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(overall.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(overall)

    pc = norm_pair(args.parent, args.child)
    pl_pc = plink.get(pc, {})
    pc_row = {
        "parent": args.parent,
        "child": args.child,
        "beagle_segments": len(beagle.get(pc, [])),
        "germline_segments": len(germ.get(pc, [])),
        "beagle_bp": total_bp(beagle.get(pc, [])),
        "germline_bp": total_bp(germ.get(pc, [])),
        "beagle_chr_prop": round(total_bp(beagle.get(pc, [])) / args.chr_length, 6),
        "germline_chr_prop": round(total_bp(germ.get(pc, [])) / args.chr_length, 6),
        "plink_pi_hat": pl_pc.get("pi_hat", ""),
        "plink_z0": pl_pc.get("z0", ""),
        "plink_z1": pl_pc.get("z1", ""),
        "plink_z2": pl_pc.get("z2", ""),
    }
    with open(outdir / "parent_child_check.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pc_row.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(pc_row)

    print("Wrote:", outdir / "overall_metrics.tsv")
    print("Wrote:", outdir / "pair_metrics.tsv")
    print("Wrote:", outdir / "parent_child_check.tsv")
    print("\n=== Summary ===")
    for k, v in overall.items():
        print(f"  {k}: {v}")
    print("\n=== Parent-Child Validation ===")
    for k, v in pc_row.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
