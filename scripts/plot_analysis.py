#!/usr/bin/env python3
"""Generate analysis figures for the GERMLINE vs Beagle IBD comparison.

Produces (in --outdir):
  genome_karyogram_{parent}_{child}.png  - per-chromosome IBD karyogram
  seg_length_hist.png                    - segment length distribution
  method_overlap.png                     - per-chromosome IBD coverage bar chart
"""
import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# hg19 chromosome lengths (bp) for all autosomes
HG19_CHR_LENGTHS = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
    6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
    11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
    16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51244237,
}

COLOR_GERM = "#b08968"
COLOR_BEAGLE = "#0b7285"


def norm_pair(a, b):
    return (a, b) if a <= b else (b, a)


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    m = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return [(a, b) for a, b in m]


def total_bp(ivs):
    return sum(e - s for s, e in merge_intervals(ivs))


def _iid(hap_id):
    if hap_id.endswith(".0") or hap_id.endswith(".1"):
        return hap_id[:-2]
    if "_" in hap_id:
        parts = hap_id.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]
    return hap_id


def parse_germline(path):
    pairs = defaultdict(list)
    seg_lengths = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")

            # germline2 format: ID1 ID2 P0 P1 cM #words #gaps
            if len(cols) >= 4 and cols[2].isdigit() and cols[3].isdigit():
                iid1 = _iid(cols[0])
                iid2 = _iid(cols[1])
                start, end = int(cols[2]), int(cols[3])
                pairs[norm_pair(iid1, iid2)].append((start, end))
                seg_lengths.append(end - start)
                continue

            # legacy GERMLINE format fallback
            p1 = cols[0].split()
            p2 = cols[1].split() if len(cols) > 1 else []
            if len(p1) < 2 or len(p2) < 2 or len(cols) < 4:
                continue
            se = cols[3].split()
            if len(se) < 2:
                continue
            iid1 = _iid(p1[1])
            iid2 = _iid(p2[1])
            start, end = int(se[0]), int(se[1])
            pairs[norm_pair(iid1, iid2)].append((start, end))
            seg_lengths.append(end - start)
    return pairs, seg_lengths


def parse_beagle(path):
    pairs = defaultdict(list)
    seg_lengths = []
    with gzip.open(path, "rb") as f:
        for line in f:
            line = line.decode().strip()
            if not line:
                continue
            sid1, _h1, sid2, _h2, _chrom, start, end, _lod = line.split("\t")
            pairs[norm_pair(sid1, sid2)].append((int(start), int(end)))
            seg_lengths.append(int(end) - int(start))
    return pairs, seg_lengths


# ── Figure 1: Genome-wide karyogram for the parent-child pair ────────────────

def plot_genome_karyogram(germ_by_chr, beagle_by_chr, chromosomes, id1, id2, outpath):
    """One row per chromosome showing GERMLINE and Beagle IBD segments."""
    pair = norm_pair(id1, id2)
    n = len(chromosomes)

    fig, axes = plt.subplots(n, 1, figsize=(12, max(4, n * 0.9)),
                             sharex=False, squeeze=False)

    for row, chrom in enumerate(chromosomes):
        ax = axes[row][0]
        chr_len = HG19_CHR_LENGTHS.get(int(chrom), 1)
        germ_segs = merge_intervals(germ_by_chr[chrom].get(pair, []))
        beagle_segs = merge_intervals(beagle_by_chr[chrom].get(pair, []))

        # Chromosome backbone
        ax.barh(0.7, chr_len, left=0, height=0.06, color="#e0e0e0", zorder=0)
        ax.barh(0.3, chr_len, left=0, height=0.06, color="#e0e0e0", zorder=0)

        for s, e in germ_segs:
            ax.barh(0.7, e - s, left=s, height=0.18,
                    color=COLOR_GERM, alpha=0.9, edgecolor="none", zorder=1)
        for s, e in beagle_segs:
            ax.barh(0.3, e - s, left=s, height=0.18,
                    color=COLOR_BEAGLE, alpha=0.9, edgecolor="none", zorder=1)

        ax.set_xlim(0, chr_len)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.3, 0.7])
        ax.set_yticklabels([
            f"Beagle ({len(beagle_segs)} segs, {total_bp(beagle_segs)/1e6:.1f} Mb)",
            f"GERMLINE ({len(germ_segs)} segs, {total_bp(germ_segs)/1e6:.1f} Mb)",
        ], fontsize=7)
        ax.set_ylabel(f"Chr{chrom}", fontsize=8, rotation=0, labelpad=30, va="center")
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(axis="x", alpha=0.15, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # x-axis label only on bottom panel
        if row < n - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Position (bp)", fontsize=8)

    germ_patch = mpatches.Patch(color=COLOR_GERM, label="GERMLINE phased")
    beagle_patch = mpatches.Patch(color=COLOR_BEAGLE, label="Beagle")
    fig.legend(handles=[germ_patch, beagle_patch], loc="upper right",
               fontsize=8, bbox_to_anchor=(0.99, 0.99))
    fig.suptitle(f"IBD Segments: {id1} vs {id2}  (ASW, 1000G Phase 3)",
                 fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── Figure 2: Segment length histogram ───────────────────────────────────────

def plot_seg_length_hist(germ_lengths, beagle_lengths, outpath):
    fig, ax = plt.subplots(figsize=(8, 4))
    bins = [b / 1e6 for b in range(0, 12_000_000, 500_000)]
    ax.hist([l / 1e6 for l in germ_lengths], bins=bins,
            alpha=0.65, color=COLOR_GERM, edgecolor="white", linewidth=0.4,
            label=f"GERMLINE phased (n={len(germ_lengths)})")
    ax.hist([l / 1e6 for l in beagle_lengths], bins=bins,
            alpha=0.65, color=COLOR_BEAGLE, edgecolor="white", linewidth=0.4,
            label=f"Beagle (n={len(beagle_lengths)})")
    ax.set_xlabel("IBD Segment Length (Mb)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of IBD Segment Lengths: GERMLINE vs Beagle")
    ax.legend()
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── Figure 3: Per-chromosome method overlap bar chart ────────────────────────

def plot_method_overlap(germ_by_chr, beagle_by_chr, chromosomes, outpath):
    def total_chr_bp(pairs_dict):
        return sum(total_bp(segs) for segs in pairs_dict.values()) / 1e6

    chrom_labels = [f"Chr{c}" for c in chromosomes]
    germ_totals = [total_chr_bp(germ_by_chr[c]) for c in chromosomes]
    beagle_totals = [total_chr_bp(beagle_by_chr[c]) for c in chromosomes]

    x = range(len(chromosomes))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, len(chromosomes) * 1.2), 4))
    bars_g = ax.bar([i - width/2 for i in x], germ_totals, width,
                    color=COLOR_GERM, alpha=0.85, label="GERMLINE phased")
    bars_b = ax.bar([i + width/2 for i in x], beagle_totals, width,
                    color=COLOR_BEAGLE, alpha=0.85, label="Beagle")

    for bar in bars_g:
        v = bar.get_height()
        if v > 0:
            ax.text(bar.get_x() + bar.get_width()/2, v + 0.5,
                    f"{v:.0f}", ha="center", va="bottom", fontsize=7)
    for bar in bars_b:
        v = bar.get_height()
        if v > 0:
            ax.text(bar.get_x() + bar.get_width()/2, v + 0.5,
                    f"{v:.0f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(list(x))
    ax.set_xticklabels(chrom_labels)
    ax.set_ylabel("Total IBD Detected (Mb)")
    ax.set_title("IBD Coverage by Method per Chromosome (ASW Population)")
    ax.legend()
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chromosomes", nargs="+", required=True)
    p.add_argument("--germline-matches", nargs="+", required=True)
    p.add_argument("--beagle-ibds", nargs="+", required=True)
    p.add_argument("--parent", default="NA20317")
    p.add_argument("--child", default="NA20318")
    p.add_argument("--outdir", default="results/figures")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    germ_by_chr = {}
    beagle_by_chr = {}
    all_germ_lengths = []
    all_beagle_lengths = []

    for chrom, gpath, bpath in zip(args.chromosomes, args.germline_matches, args.beagle_ibds):
        print(f"Parsing chr{chrom}...")
        germ_pairs, germ_lens = parse_germline(gpath)
        beagle_pairs, beagle_lens = parse_beagle(bpath)
        germ_by_chr[chrom] = germ_pairs
        beagle_by_chr[chrom] = beagle_pairs
        all_germ_lengths.extend(germ_lens)
        all_beagle_lengths.extend(beagle_lens)

    plot_genome_karyogram(
        germ_by_chr, beagle_by_chr, args.chromosomes,
        args.parent, args.child,
        outdir / f"genome_karyogram_{args.parent}_{args.child}.png",
    )
    plot_seg_length_hist(all_germ_lengths, all_beagle_lengths,
                         outdir / "seg_length_hist.png")
    plot_method_overlap(germ_by_chr, beagle_by_chr, args.chromosomes,
                        outdir / "method_overlap.png")


if __name__ == "__main__":
    main()
