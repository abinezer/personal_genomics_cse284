#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import matplotlib.pyplot as plt


def norm_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def read_beagle(path: Path, target: tuple[str, str]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            sid1, _h1, sid2, _h2, _chrom, start, end, _lod = line.strip().split("\t")
            if norm_pair(sid1, sid2) == target:
                out.append((int(start), int(end)))
    return out


def _strip_hap(hap_id: str) -> str:
    """Strip '_1' / '_2' haplotype suffix from GERMLINE phased IDs."""
    if "_" in hap_id:
        parts = hap_id.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]
    return hap_id


def read_germline(path: Path, target: tuple[str, str]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split("\t")
            sid1 = _strip_hap(cols[0].split()[1])
            sid2 = _strip_hap(cols[1].split()[1])
            if norm_pair(sid1, sid2) == target:
                s, e = cols[3].split()
                out.append((int(s), int(e)))
    return out


def draw_segments(ax, segs: list[tuple[int, int]], y: float, color: str, label: str):
    for s, e in segs:
        ax.plot([s, e], [y, y], color=color, linewidth=8, solid_capstyle="butt")
    ax.text(0, y + 0.05, f"{label} (n={len(segs)})", fontsize=9, ha="left", va="bottom")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--beagle-ibd", required=True)
    p.add_argument("--germline-match", required=True)
    p.add_argument("--id1", required=True)
    p.add_argument("--id2", required=True)
    p.add_argument("--chr-length", type=int, default=51_244_237)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    target = norm_pair(args.id1, args.id2)
    beagle = read_beagle(Path(args.beagle_ibd), target)
    germ = read_germline(Path(args.germline_match), target)

    fig, ax = plt.subplots(figsize=(10, 2.8))
    draw_segments(ax, beagle, 1.0, "#0b7285", "Beagle")
    draw_segments(ax, germ, 0.4, "#b08968", "GERMLINE")
    ax.set_xlim(0, args.chr_length)
    ax.set_ylim(0, 1.4)
    ax.set_yticks([])
    ax.set_xlabel("Chr22 Position (bp)")
    ax.set_title(f"IBD Segments for {target[0]} vs {target[1]}")
    ax.grid(axis="x", alpha=0.2, linewidth=0.6)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(args.out)


if __name__ == "__main__":
    main()
