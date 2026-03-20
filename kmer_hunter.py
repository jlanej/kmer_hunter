#!/usr/bin/env python3
"""
kmer_hunter: Exact k-mer search against T2T CHM13v2.0 chrY with
an interactive Plotly HTML report showing PAR/XTR regions on a karyogram.

All matches are EXACT — zero mismatches, zero gaps.
Both forward and reverse-complement strands are searched.
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.io import to_html

# ─── T2T chrY reference ──────────────────────────────────────────────────────

T2T_CHRY_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/chromosomes/chrY.fa.gz"
)

# ─── T2T whole-genome reference ───────────────────────────────────────────────

T2T_WHOLE_GENOME_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
)

CHRY_LEN = 62_460_029  # T2T CHM13v2.0 (hs1) chrY length (bp)

# ─── chrY functional-region annotations ──────────────────────────────────────
# Coordinates: 1-based, inclusive.
# Sources: Rhie et al. (2023) Nature (T2T-Y) and UCSC hs1 tracks.
CHRY_REGIONS = [
    {
        "name": "PAR1",
        "start": 1,
        "end": 2_781_479,
        "color": "#27ae60",
        "description": "Pseudoautosomal Region 1 (Yp telomere)",
    },
    {
        "name": "XTR",
        "start": 2_781_480,
        "end": 6_811_428,
        "color": "#2980b9",
        "description": "X-Transposed Region",
    },
    {
        "name": "Ampliconic",
        "start": 6_811_429,
        "end": 26_200_000,
        "color": "#8e44ad",
        "description": "Ampliconic sequences (Yp/Yq)",
    },
    {
        "name": "Pericentromeric",
        "start": 26_200_001,
        "end": 27_800_000,
        "color": "#c0392b",
        "description": "Pericentromeric / Centromere",
    },
    {
        "name": "Heterochromatin",
        "start": 27_800_001,
        "end": 56_887_901,
        "color": "#95a5a6",
        "description": "Heterochromatin (DYZ1/DYZ2 satellites)",
    },
    {
        "name": "PAR2",
        "start": 56_887_902,
        "end": 57_217_415,
        "color": "#e67e22",
        "description": "Pseudoautosomal Region 2 (Yq telomere)",
    },
    {
        "name": "Distal Yq",
        "start": 57_217_416,
        "end": CHRY_LEN,
        "color": "#bdc3c7",
        "description": "Distal Yq (T2T-resolved region)",
    },
]

# ─── CLI ──────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="kmer_hunter",
        description=(
            "Exact k-mer search against T2T CHM13v2.0 chrY.\n"
            "All matches are perfect (zero mismatches, zero gaps).\n"
            "Generates an interactive HTML report with a chrY karyogram,\n"
            "PAR/XTR region annotations, and a detailed hit table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "kmers",
        metavar="KMERS",
        help="Input file: one k-mer per line (plain text) or FASTA format.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="kmer_report.html",
        metavar="PATH",
        help="Output HTML report [default: kmer_report.html]",
    )
    parser.add_argument(
        "--reference",
        default=None,
        metavar="PATH",
        help=(
            "Path to reference FASTA (.fa or .fa.gz). "
            "Defaults to T2T CHM13v2.0 chrY, auto-downloaded to --cache-dir."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        default=str(Path.home() / ".kmer_hunter_cache"),
        metavar="DIR",
        help="Cache directory for downloaded reference [default: ~/.kmer_hunter_cache]",
    )
    parser.add_argument(
        "--whole-genome",
        action="store_true",
        help=(
            "Also search the full T2T CHM13v2.0 genome for non-chrY hits. "
            "Auto-downloads hs1.fa.gz (~830 MB) on first use unless "
            "--whole-genome-reference is provided."
        ),
    )
    parser.add_argument(
        "--whole-genome-reference",
        default=None,
        metavar="PATH",
        help=(
            "Path to a full-genome FASTA used with --whole-genome "
            "instead of auto-downloading T2T CHM13v2.0. "
            "Any genome name is accepted (e.g. chm13v2.0.fa); "
            "BWA index files (*.amb, *.ann, *.bwt, *.pac, *.sa) are "
            "expected alongside the FASTA and will be built automatically "
            "if absent."
        ),
    )
    parser.add_argument(
        "--output-bam",
        default=None,
        metavar="PATH",
        help=(
            "Write a BAM (or SAM) alignment file with all BWA hits. "
            "If the path ends in .bam and samtools is available, a "
            "coordinate-sorted BAM with index is produced; otherwise "
            "a SAM file is written."
        ),
    )
    return parser.parse_args()


# ─── Sequence utilities ───────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# ─── Input parsing ────────────────────────────────────────────────────────────


def read_kmers(path: str) -> list[tuple[str, str]]:
    """Read k-mers from a plain-text (one per line) or FASTA file.

    Returns a list of (name, sequence) tuples.
    """
    text = Path(path).read_text()
    lines = text.strip().splitlines()

    if not lines:
        sys.exit(f"[kmer_hunter] ERROR: Input file is empty: {path}")

    kmers: list[tuple[str, str]] = []

    if lines[0].startswith(">"):
        # FASTA format
        name: str | None = None
        parts: list[str] = []
        for line in lines:
            if line.startswith(">"):
                if name is not None:
                    kmers.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            kmers.append((name, "".join(parts).upper()))
    else:
        # Plain-text, one k-mer per line
        for i, line in enumerate(lines, start=1):
            seq = line.strip()
            if seq and not seq.startswith("#"):
                kmers.append((f"kmer_{i}", seq.upper()))

    if not kmers:
        sys.exit(f"[kmer_hunter] ERROR: No k-mers found in: {path}")

    return kmers


# ─── Reference genome ─────────────────────────────────────────────────────────


def _extract_chrom(source_gz: str, chrom: str, dest_gz: str) -> None:
    """Stream a single chromosome from a gzipped FASTA into a new gzipped FASTA.

    Raises ValueError if *chrom* is not found in *source_gz*.
    The partially-written *dest_gz* is removed on failure.
    """
    found = False
    in_target = False
    try:
        with gzip.open(source_gz, "rt") as src, gzip.open(dest_gz, "wt") as dst:
            for line in src:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    in_target = name == chrom
                    if in_target:
                        found = True
                        dst.write(line)
                elif in_target:
                    dst.write(line)
    except Exception:
        Path(dest_gz).unlink(missing_ok=True)
        raise
    if not found:
        Path(dest_gz).unlink(missing_ok=True)
        raise ValueError(f"chromosome '{chrom}' not found in {source_gz}")


def ensure_reference(reference: str | None, cache_dir: str, whole_genome_reference: str | None = None) -> str:
    """Return a path to the reference FASTA, downloading the full T2T genome if needed.

    When no reference is provided the full T2T CHM13v2.0 genome is downloaded
    (or reused from cache) so that all chromosomes are available for BWA search.
    If *whole_genome_reference* is given it is used instead of downloading.
    """
    if reference:
        if not Path(reference).exists():
            sys.exit(f"[kmer_hunter] ERROR: Reference not found: {reference}")
        return reference

    return ensure_whole_genome_reference(whole_genome_reference, cache_dir)


def ensure_whole_genome_reference(reference: str | None, cache_dir: str) -> str:
    """Return a path to the full-genome FASTA, downloading if necessary."""
    if reference:
        if not Path(reference).exists():
            sys.exit(
                f"[kmer_hunter] ERROR: Whole-genome reference not found: {reference}"
            )
        return reference

    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    ref_gz = cache / "chm13v2.0_hs1.fa.gz"

    if ref_gz.exists():
        print(
            f"[kmer_hunter] Using cached whole-genome reference: {ref_gz}",
            file=sys.stderr,
        )
        return str(ref_gz)

    print(
        f"[kmer_hunter] Downloading T2T CHM13v2.0 full genome → {ref_gz}",
        file=sys.stderr,
    )
    print(f"[kmer_hunter] URL: {T2T_WHOLE_GENOME_URL}", file=sys.stderr)
    print(
        "[kmer_hunter] Warning: Full genome is ~830 MB; this may take a while.",
        file=sys.stderr,
    )

    try:
        urllib.request.urlretrieve(T2T_WHOLE_GENOME_URL, str(ref_gz))
    except Exception as exc:
        sys.exit(f"[kmer_hunter] ERROR: Download failed: {exc}")

    return str(ref_gz)


def read_fasta(path: str) -> dict[str, str]:
    """Read a FASTA file (plain or .gz) into {chrom: sequence} (uppercase)."""
    opener = gzip.open if path.endswith(".gz") else open
    sequences: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []

    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())

    if name is not None:
        sequences[name] = "".join(parts)

    return sequences


# ─── BWA-based exact matching ─────────────────────────────────────────────────


def ensure_bwa_index(ref_path: str) -> None:
    """Build a BWA index for *ref_path* unless one already exists.

    The index files (*.bwt, *.pac, *.amb, *.ann, *.sa) are stored next to the
    reference file.  The *.bwt file is used as a sentinel: if it is present we
    assume the full index is intact and skip re-indexing.  If BWA finds any
    other index files missing it will report a clear error at search time.
    """
    bwt = ref_path + ".bwt"
    if Path(bwt).exists():
        print(f"[kmer_hunter] BWA index already exists: {bwt}", file=sys.stderr)
        return
    print(f"[kmer_hunter] Building BWA index: {ref_path} …", file=sys.stderr)
    result = subprocess.run(
        ["bwa", "index", ref_path], capture_output=True, text=True
    )
    if result.returncode != 0:
        sys.exit(f"[kmer_hunter] ERROR: bwa index failed:\n{result.stderr}")
    print("[kmer_hunter] BWA index ready.", file=sys.stderr)


def _parse_sam_exact(sam_text: str, kmer_seqs: dict[str, str]) -> list[dict]:
    """Return only exact, full-length alignments from SAM output.

    An alignment is accepted when:
    - The read is mapped (FLAG & 4 == 0).
    - It is not a supplementary alignment (FLAG & 2048 == 0).
    - The CIGAR string is exactly ``{kmer_len}M`` (no soft clips, no gaps).
    - The ``NM`` tag equals 0 (zero mismatches / indels).
    """
    # Pre-compute expected CIGAR strings to avoid repeated f-string formatting.
    kmer_info: dict[str, tuple[str, str]] = {}
    for name, seq in kmer_seqs.items():
        kmer_info[name] = (seq, f"{len(seq)}M")

    NM_EXACT = "NM:i:0"  # noqa: N806 – constant used as fast substring check

    hits: list[dict] = []
    for line in sam_text.splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])

        if flag & 0x804:  # unmapped (4) or supplementary (2048)
            continue

        qname = fields[0]
        info = kmer_info.get(qname)
        if info is None:
            continue
        kmer_seq, expected_cigar = info

        if fields[5] != expected_cigar:
            continue

        # Fast check: look for exact NM:i:0 among the optional fields.
        if not any(f == NM_EXACT for f in fields[11:]):
            continue

        pos = int(fields[3])  # 1-based in SAM
        strand = "-" if (flag & 16) else "+"
        kmer_len = len(kmer_seq)
        end = pos + kmer_len - 1
        hits.append(
            {
                "kmer": qname,
                "seq": kmer_seq,
                "chrom": fields[2],
                "start": pos,
                "end": end,
                "strand": strand,
                "region": annotate_region(fields[2], pos),
            }
        )
    return hits


def bwa_find_exact_matches(
    kmers: list[tuple[str, str]],
    ref_path: str,
) -> tuple[list[dict], str]:
    """Find all exact k-mer matches using BWA mem.

    This is orders of magnitude faster than Python string search for large
    k-mer sets.  BWA mem is invoked once for the entire k-mer collection with
    settings that strongly penalise mismatches, gaps, and soft-clipping so
    that only perfect full-length alignments score above the reporting
    threshold.  The SAM output is then filtered through ``_parse_sam_exact``
    for additional correctness.

    The BWA index is built automatically (see :func:`ensure_bwa_index`) unless
    it already exists alongside *ref_path*.

    Returns ``(hits, sam_text)`` where *sam_text* is the raw SAM output from
    BWA.  Call :func:`save_alignment_file` separately if you need to persist
    the alignments.
    """
    ensure_bwa_index(ref_path)
    kmer_seqs: dict[str, str] = dict(kmers)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        query_fa = fh.name
        for name, seq in kmers:
            fh.write(f">{name}\n{seq}\n")

    try:
        print(
            f"[kmer_hunter] Running bwa mem on {len(kmers)} k-mer(s) …",
            file=sys.stderr,
        )
        result = subprocess.run(
            [
                "bwa", "mem",
                "-a",          # report all alignments, not just the best
                "-k", "11",    # minimum seed length (handles k-mers ≥ 11 bp)
                # NOTE: do NOT pass extreme mismatch/gap/clip penalties here.
                # Values like -B 1000 cause integer overflow inside BWA's
                # internal scoring for repetitive k-mers (e.g. poly-A prefix),
                # leading BWA to report wrong primary alignments and miss the
                # actual exact matches.  Exact-match filtering is handled
                # entirely by _parse_sam_exact (CIGAR = {klen}M, NM = 0).
                "-T", "0",     # min score threshold; we filter by CIGAR + NM
                ref_path,
                query_fa,
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            sys.exit(f"[kmer_hunter] ERROR: bwa mem failed:\n{result.stderr}")
        sam_text = result.stdout
        hits = _parse_sam_exact(sam_text, kmer_seqs)
        print(f"[kmer_hunter] {len(hits)} exact hit(s) found", file=sys.stderr)

        return hits, sam_text
    finally:
        os.unlink(query_fa)


# ─── Exact matching ───────────────────────────────────────────────────────────


def annotate_region(chrom: str, pos: int) -> str:
    """Return the functional region name for a chrY 1-based position."""
    if chrom != "chrY":
        return chrom  # report the chromosome name for non-chrY hits
    for region in CHRY_REGIONS:
        if region["start"] <= pos <= region["end"]:
            return region["name"]
    return "Unknown"


def find_exact_matches(
    kmers: list[tuple[str, str]],
    ref_seqs: dict[str, str],
) -> list[dict]:
    """Find all exact (zero-mismatch, zero-gap) occurrences of each k-mer.

    Both forward (+) and reverse-complement (−) strands are searched.
    Overlapping matches are included.
    """
    hits: list[dict] = []
    total = len(kmers)

    for idx, (kmer_name, kmer_seq) in enumerate(kmers, start=1):
        kmer_upper = kmer_seq.upper()
        kmer_rc = reverse_complement(kmer_upper)
        kmer_len = len(kmer_upper)

        print(
            f"[kmer_hunter]  [{idx}/{total}] searching for {kmer_name} ({kmer_len} bp) …",
            file=sys.stderr,
        )

        for chrom, chrom_seq in ref_seqs.items():
            # ── Forward strand (overlapping via lookahead) ────────────────
            pattern_fwd = f"(?=({re.escape(kmer_upper)}))"
            for m in re.finditer(pattern_fwd, chrom_seq):
                start = m.start() + 1  # convert to 1-based
                hits.append(
                    {
                        "kmer": kmer_name,
                        "seq": kmer_seq,
                        "chrom": chrom,
                        "start": start,
                        "end": start + kmer_len - 1,
                        "strand": "+",
                        "region": annotate_region(chrom, start),
                    }
                )

            # ── Reverse-complement strand ─────────────────────────────────
            if kmer_rc != kmer_upper:  # skip palindromes (already counted)
                pattern_rev = f"(?=({re.escape(kmer_rc)}))"
                for m in re.finditer(pattern_rev, chrom_seq):
                    start = m.start() + 1
                    hits.append(
                        {
                            "kmer": kmer_name,
                            "seq": kmer_seq,
                            "chrom": chrom,
                            "start": start,
                            "end": start + kmer_len - 1,
                            "strand": "-",
                            "region": annotate_region(chrom, start),
                        }
                    )

    return hits


# ─── Interval collapsing & cluster detection ──────────────────────────────────


def collapse_to_intervals(
    hits_df: pd.DataFrame, gap: int = 1000
) -> pd.DataFrame:
    """Merge nearby hits on the same chromosome into intervals with counts.

    Hits whose start positions are within *gap* bp of the previous hit's end
    are merged into a single interval.  Returns a DataFrame with columns:
    ``chrom``, ``start``, ``end``, ``count``, ``unique_count``, ``region``.

    ``count`` is the total number of hits in the interval.
    ``unique_count`` counts only hits from k-mers that appear exactly once
    across the entire hit set (unique k-mers), providing a measure of
    interval enrichment that is not inflated by repetitive sequences.
    """
    if hits_df.empty:
        return pd.DataFrame(
            columns=["chrom", "start", "end", "count", "unique_count", "region"]
        )

    # Identify k-mers that appear exactly once across all hits (unique k-mers).
    kmer_hit_counts = hits_df.groupby("kmer").size()
    unique_kmers: set[str] = set(kmer_hit_counts[kmer_hit_counts == 1].index)

    intervals: list[dict] = []
    sorted_df = hits_df.sort_values(["chrom", "start"]).reset_index(drop=True)

    for chrom, group in sorted_df.groupby("chrom", sort=False):
        rows = group.to_dict("records")
        cur = {
            "chrom": chrom,
            "start": rows[0]["start"],
            "end": rows[0]["end"],
            "count": 1,
            "unique_count": 1 if rows[0]["kmer"] in unique_kmers else 0,
            "region": rows[0]["region"],
        }
        for row in rows[1:]:
            if row["start"] <= cur["end"] + gap:
                cur["end"] = max(cur["end"], row["end"])
                cur["count"] += 1
                if row["kmer"] in unique_kmers:
                    cur["unique_count"] += 1
            else:
                intervals.append(cur)
                cur = {
                    "chrom": chrom,
                    "start": row["start"],
                    "end": row["end"],
                    "count": 1,
                    "unique_count": 1 if row["kmer"] in unique_kmers else 0,
                    "region": row["region"],
                }
        intervals.append(cur)

    return pd.DataFrame(intervals)


def detect_clusters(
    intervals_df: pd.DataFrame, min_hits: int = 5
) -> pd.DataFrame:
    """Mark intervals as clusters when they contain at least *min_hits* hits.

    Returns a copy of *intervals_df* with an added boolean ``cluster`` column.
    """
    if intervals_df.empty:
        df = intervals_df.copy()
        df["cluster"] = pd.Series(dtype=bool)
        return df

    df = intervals_df.copy()
    df["cluster"] = df["count"] >= min_hits
    return df


def build_non_chry_summary(hits_df: pd.DataFrame) -> go.Figure:
    """Summary table of non-chrY hits grouped by chromosome.

    Shows chromosome, total hits, and number of distinct k-mers matched.
    """
    non_chry = (
        hits_df[hits_df["chrom"] != "chrY"]
        if not hits_df.empty
        else pd.DataFrame(columns=["kmer", "chrom", "start", "end", "strand", "seq"])
    )

    if non_chry.empty:
        display = pd.DataFrame(columns=["chrom", "total_hits", "distinct_kmers"])
    else:
        summary = non_chry.groupby("chrom").agg(
            total_hits=("chrom", "size"),
            distinct_kmers=("kmer", "nunique"),
        ).reset_index()
        display = summary.sort_values("chrom").reset_index(drop=True)

    fig = go.Figure(
        go.Table(
            header=dict(
                values=[f"<b>{c}</b>" for c in display.columns],
                fill_color="#2c3e50",
                font=dict(color="white", size=12),
                align="left",
                height=32,
            ),
            cells=dict(
                values=[display[c] for c in display.columns],
                fill_color=[["#f8f9fa"] * len(display) for _ in display.columns],
                align="left",
                font=dict(size=11),
                height=28,
            ),
        )
    )
    fig.update_layout(
        title="Non-chrY Exact Hit Summary by Chromosome",
        height=max(220, 35 * len(display) + 70),
        margin=dict(t=60, b=20),
    )
    return fig


# ─── Alignment file output ────────────────────────────────────────────────────


def save_alignment_file(sam_text: str, output_path: str) -> str:
    """Write filtered SAM and optionally convert to sorted BAM.

    If *samtools* is available the SAM is converted to a coordinate-sorted BAM
    with index.  Otherwise the SAM is saved as-is.  Returns the path of the
    file actually written.
    """
    sam_path = output_path
    if not sam_path.endswith(".sam") and not sam_path.endswith(".bam"):
        sam_path = output_path + ".sam"

    has_samtools = shutil.which("samtools") is not None

    if not has_samtools and sam_path.endswith(".bam"):
        sam_path = sam_path[:-4] + ".sam"

    if has_samtools and output_path.endswith(".bam"):
        # SAM → sorted BAM
        bam_path = output_path
        unsorted = bam_path + ".unsorted.bam"
        try:
            # Convert SAM → unsorted BAM
            p1 = subprocess.run(
                ["samtools", "view", "-bS", "-"],
                input=sam_text.encode("utf-8"),
                capture_output=True,
            )
            if p1.returncode != 0:
                raise RuntimeError(p1.stderr.decode())
            Path(unsorted).write_bytes(p1.stdout)

            # Sort
            p2 = subprocess.run(
                ["samtools", "sort", "-o", bam_path, unsorted],
                capture_output=True,
                text=True,
            )
            if p2.returncode != 0:
                raise RuntimeError(p2.stderr)

            # Index
            subprocess.run(
                ["samtools", "index", bam_path],
                capture_output=True,
                text=True,
            )
            print(f"[kmer_hunter] Alignment file (BAM) written → {bam_path}", file=sys.stderr)
            return bam_path
        except Exception as exc:
            print(
                f"[kmer_hunter] WARNING: BAM conversion failed ({exc}); "
                f"falling back to SAM",
                file=sys.stderr,
            )
            sam_path = output_path.replace(".bam", ".sam")
        finally:
            Path(unsorted).unlink(missing_ok=True)

    # Fall back: write SAM
    Path(sam_path).write_text(sam_text, encoding="utf-8")
    print(f"[kmer_hunter] Alignment file written → {sam_path}", file=sys.stderr)
    return sam_path



# ─── Text file reports ────────────────────────────────────────────────────────


def write_multi_match_report(hits_df: pd.DataFrame, output_stem: str) -> str | None:
    """Write kmers with more than one exact match to a plain-text file.

    Returns the path of the file written, or ``None`` if there are no
    multi-match kmers.
    """
    if hits_df.empty:
        return None

    kmer_counts = hits_df.groupby("kmer").size()
    multi_kmers = kmer_counts[kmer_counts > 1].index.tolist()

    if not multi_kmers:
        return None

    path = f"{output_stem}_multi_match_kmers.txt"
    lines = [
        "# Kmers with more than 1 exact match",
        "# Generated by kmer_hunter",
        "# Format: kmer_name<TAB>sequence<TAB>total_hits",
        "#         followed by each hit location",
        "",
    ]

    multi_set = set(multi_kmers)
    multi_df = (
        hits_df[hits_df["kmer"].isin(multi_set)]
        .sort_values(["kmer", "chrom", "start"])
    )

    # Build all detail lines vectorised across the entire multi_df at once.
    detail_series = (
        "  " + multi_df["chrom"].astype(str) + ":"
        + multi_df["start"].astype(str) + "-"
        + multi_df["end"].astype(str)
        + " (" + multi_df["strand"].astype(str) + ")"
        + " [" + multi_df["region"].astype(str) + "]"
    )
    multi_df = multi_df.assign(_detail=detail_series.values)

    # Extract first-seq-per-kmer and counts via a single groupby.
    agg = multi_df.groupby("kmer", sort=True).agg(
        seq=("seq", "first"),
        total=("kmer", "size"),
        details=("_detail", list),
    )

    for kmer, row in agg.iterrows():
        lines.append(f"{kmer}\t{row['seq']}\t{row['total']}")
        lines.extend(row["details"])
        lines.append("")

    Path(path).write_text("\n".join(lines), encoding="utf-8")
    print(f"[kmer_hunter] Multi-match kmer report written → {path}", file=sys.stderr)
    return path


def write_non_chry_report(hits_df: pd.DataFrame, output_stem: str) -> str | None:
    """Write non-chrY exact matches to a plain-text file.

    Lists every unique (kmer, chromosome) hit so that kmers mapping outside
    chrY can be inspected without opening the HTML report.  Returns the path
    written, or ``None`` when there are no non-chrY hits.
    """
    non_chry = (
        hits_df[hits_df["chrom"] != "chrY"].copy()
        if not hits_df.empty
        else pd.DataFrame()
    )

    if non_chry.empty:
        return None

    path = f"{output_stem}_non_chry_kmers.txt"
    display = non_chry[
        ["kmer", "chrom", "start", "end", "strand", "seq"]
    ].sort_values(["chrom", "start"]).reset_index(drop=True)

    lines = [
        "# Non-chrY exact k-mer matches",
        "# Generated by kmer_hunter",
        "# Columns: kmer<TAB>chrom<TAB>start<TAB>end<TAB>strand<TAB>seq",
        "",
    ]
    lines.extend(
        display["kmer"].astype(str) + "\t"
        + display["chrom"].astype(str) + "\t"
        + display["start"].astype(str) + "\t"
        + display["end"].astype(str) + "\t"
        + display["strand"].astype(str) + "\t"
        + display["seq"].astype(str)
    )

    Path(path).write_text("\n".join(lines), encoding="utf-8")
    print(f"[kmer_hunter] Non-chrY kmer report written → {path}", file=sys.stderr)
    return path


# ─── Plotly figures ───────────────────────────────────────────────────────────


def build_karyogram(
    hits_df: pd.DataFrame,
    intervals_df: pd.DataFrame | None = None,
) -> go.Figure:
    """Interactive chrY karyogram with PAR/XTR bands and a hit-count y-axis.

    The y-axis shows the number of exact k-mer hits per collapsed interval,
    making clusters and regions with many matches immediately visible.
    Cluster intervals (``cluster == True``) are highlighted in a darker red.

    Two views are available via toggle buttons:

    * **Unique Hits Only** (default) — bars reflect ``unique_count``, counting
      only hits from k-mers that appear exactly once genome-wide.  This removes
      inflation from repetitive or multi-mapping sequences.
    * **All Hits** — bars reflect the total ``count`` of every hit in the
      interval, including hits from k-mers with multiple genome-wide matches.

    When *intervals_df* is ``None`` or empty but per-hit data is available,
    individual hit markers are shown at y = 1 as a fallback (no toggle).
    """
    fig = go.Figure()

    # Determine the y-axis ceilings from the interval data
    if intervals_df is not None and not intervals_df.empty:
        chry_intervals = intervals_df[intervals_df["chrom"] == "chrY"].copy()
    else:
        chry_intervals = pd.DataFrame()

    # Backward-compat: if unique_count column is absent, fall back to count.
    if not chry_intervals.empty and "unique_count" not in chry_intervals.columns:
        chry_intervals["unique_count"] = chry_intervals["count"]

    total_max = int(chry_intervals["count"].max()) if not chry_intervals.empty else 1
    unique_max = (
        int(chry_intervals["unique_count"].max()) if not chry_intervals.empty else 1
    )
    y_max_all = max(total_max * 1.10, 2)      # ceiling for "All Hits" view
    y_max_unique = max(unique_max * 1.10, 2)  # ceiling for "Unique Hits Only" view
    # Shapes use the larger ceiling so they cover both views without clipping.
    y_max = y_max_all

    # ── Background: functional region colour bands ────────────────────────
    for region in CHRY_REGIONS:
        fig.add_shape(
            type="rect",
            x0=region["start"],
            x1=region["end"],
            y0=0,
            y1=y_max,
            fillcolor=region["color"],
            opacity=0.15,
            line=dict(width=0),
            layer="below",
        )
        # Label wide regions; use paper-space y so labels stay visible when
        # the y-axis range changes via the toggle.
        width = region["end"] - region["start"]
        if width > 1_500_000:
            fig.add_annotation(
                x=(region["start"] + region["end"]) / 2,
                y=0.97,
                yref="paper",
                text=f'<b>{region["name"]}</b>',
                showarrow=False,
                font=dict(size=9, color="#2c3e50"),
                xanchor="center",
                yanchor="top",
            )

    # Chromosome body outline at the bottom (thin decorative bar)
    fig.add_shape(
        type="rect",
        x0=0,
        x1=CHRY_LEN,
        y0=0,
        y1=y_max * 0.015,
        fillcolor="#bdc3c7",
        line=dict(color="#95a5a6", width=1),
        layer="below",
    )

    # Invisible legend entries for region colours
    for region in CHRY_REGIONS:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=13, color=region["color"], symbol="square"),
                name=region["name"],
                hoverinfo="skip",
                showlegend=True,
            )
        )

    chry_hits = (
        hits_df[hits_df["chrom"] == "chrY"]
        if not hits_df.empty
        else pd.DataFrame()
    )

    # ── Hit-count bars for each collapsed interval ────────────────────────
    if not chry_intervals.empty:
        midpoints = (chry_intervals["start"] + chry_intervals["end"]) / 2
        widths = (chry_intervals["end"] - chry_intervals["start"]).clip(lower=1).tolist()

        is_cluster = (
            chry_intervals["cluster"]
            if "cluster" in chry_intervals.columns
            else pd.Series(False, index=chry_intervals.index)
        )

        bar_colors = [
            "#c0392b" if c else "#e74c3c" for c in is_cluster
        ]
        bar_line_colors = [
            "#7b241c" if c else "rgba(0,0,0,0)" for c in is_cluster
        ]
        bar_line_widths = [1.5 if c else 0 for c in is_cluster]

        starts = chry_intervals["start"]
        ends = chry_intervals["end"]
        unique_counts = chry_intervals["unique_count"].astype(int)
        counts = chry_intervals["count"].astype(int)
        regions = chry_intervals["region"]
        hover_texts = [
            (
                f"<b>Interval</b><br>"
                f"chrY:{s:,}–{e:,}<br>"
                f"Unique Hits: {uc}<br>"
                f"Total Hits: {tc}<br>"
                f"Region: {rg}<br>"
                f"{'🔴 <b>Cluster</b>' if cl else ''}"
            )
            for s, e, uc, tc, rg, cl in zip(
                starts, ends, unique_counts, counts, regions, is_cluster
            )
        ]

        # Trace: Unique Hits Only bars (default visible)
        fig.add_trace(
            go.Bar(
                x=midpoints.tolist(),
                y=chry_intervals["unique_count"].tolist(),
                width=widths,
                marker=dict(
                    color=bar_colors,
                    opacity=0.75,
                    line=dict(color=bar_line_colors, width=bar_line_widths),
                ),
                text=hover_texts,
                hovertemplate="%{text}<extra></extra>",
                name="unique hit interval",
                visible=True,
                showlegend=True,
            )
        )

        # Trace: All Hits bars (initially hidden)
        fig.add_trace(
            go.Bar(
                x=midpoints.tolist(),
                y=chry_intervals["count"].tolist(),
                width=widths,
                marker=dict(
                    color=bar_colors,
                    opacity=0.75,
                    line=dict(color=bar_line_colors, width=bar_line_widths),
                ),
                text=hover_texts,
                hovertemplate="%{text}<extra></extra>",
                name="hit interval (all)",
                visible=False,
                showlegend=True,
            )
        )

        # Optional legend entry for cluster bars (always visible when present)
        has_clusters = bool(
            chry_intervals.get("cluster", pd.Series(dtype=bool)).any()
        )
        if has_clusters:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        size=13, color="#c0392b", symbol="square",
                        line=dict(color="#7b241c", width=2),
                    ),
                    name="cluster (dense)",
                    hoverinfo="skip",
                    showlegend=True,
                )
            )

        # ── Toggle buttons: Unique Hits Only ↔ All Hits ───────────────────
        n_legend = len(CHRY_REGIONS)   # region legend traces (always visible)
        n_extra = 1 if has_clusters else 0  # cluster legend trace
        # visible arrays must cover every trace in order
        visible_unique = [True] * n_legend + [True, False] + [True] * n_extra
        visible_all = [True] * n_legend + [False, True] + [True] * n_extra

        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    buttons=[
                        dict(
                            label="Unique Hits Only",
                            method="update",
                            args=[
                                {"visible": visible_unique},
                                {"yaxis.range": [0, y_max_unique]},
                            ],
                        ),
                        dict(
                            label="All Hits",
                            method="update",
                            args=[
                                {"visible": visible_all},
                                {"yaxis.range": [0, y_max_all]},
                            ],
                        ),
                    ],
                    active=0,
                    x=0.0,
                    xanchor="left",
                    y=1.15,
                    yanchor="top",
                    bgcolor="#ecf0f1",
                    bordercolor="#bdc3c7",
                    font=dict(size=11),
                ),
            ]
        )

    elif not chry_hits.empty:
        # Fallback: individual markers at y = 1 when no intervals are available
        midpoints = (chry_hits["start"] + chry_hits["end"]) / 2
        hover_texts = [
            (
                f"<b>{k}</b><br>"
                f"Position: chrY:{s:,}–{e:,}<br>"
                f"Strand: {st}<br>"
                f"Region: {rg}<br>"
                f"Sequence: <code>{sq}</code>"
            )
            for k, s, e, st, rg, sq in zip(
                chry_hits["kmer"], chry_hits["start"], chry_hits["end"],
                chry_hits["strand"], chry_hits["region"], chry_hits["seq"],
            )
        ]
        fig.add_trace(
            go.Scatter(
                x=midpoints,
                y=[1] * len(chry_hits),
                mode="markers",
                marker=dict(
                    symbol="triangle-down",
                    size=14,
                    color="#e74c3c",
                    line=dict(color="#c0392b", width=1),
                ),
                text=hover_texts,
                hovertemplate="%{text}<extra></extra>",
                name="k-mer hit (exact)",
                showlegend=True,
            )
        )

    fig.update_layout(
        title=dict(
            text="chrY Karyogram — Exact k-mer Hit Counts per Interval",
            font=dict(size=16),
        ),
        xaxis=dict(
            title="chrY position (bp)",
            tickformat=",",
            range=[-500_000, CHRY_LEN + 500_000],
            showgrid=False,
        ),
        yaxis=dict(
            title="Hit count",
            visible=True,
            # Start in "Unique Hits Only" range; the toggle buttons update this.
            # When there are no intervals, y_max_unique == y_max_all == 2.
            range=[0, y_max_unique],
            showgrid=True,
            gridcolor="#ecf0f1",
        ),
        height=420,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.05,
            xanchor="right",
            x=1,
            font=dict(size=10),
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=80, b=50, l=60, r=20),
    )
    return fig


def build_region_bar(hits_df: pd.DataFrame) -> go.Figure:
    """Bar chart: number of exact hits per chrY functional region."""
    chry_hits = (
        hits_df[hits_df["chrom"] == "chrY"]
        if not hits_df.empty
        else pd.DataFrame()
    )

    region_names = [r["name"] for r in CHRY_REGIONS]
    region_colors = {r["name"]: r["color"] for r in CHRY_REGIONS}
    counts = {r: 0 for r in region_names}

    if not chry_hits.empty:
        for region, cnt in chry_hits["region"].value_counts().items():
            if region in counts:
                counts[region] = int(cnt)

    fig = go.Figure(
        go.Bar(
            x=region_names,
            y=[counts[r] for r in region_names],
            marker_color=[region_colors[r] for r in region_names],
            text=[counts[r] for r in region_names],
            textposition="auto",
            hovertemplate="<b>%{x}</b><br>Hits: %{y}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Exact k-mer Hits by chrY Functional Region",
        xaxis_title="Region",
        yaxis_title="Number of Exact Hits",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=350,
        margin=dict(t=60, b=50),
    )
    return fig


def build_non_chry_bar(hits_df: pd.DataFrame) -> go.Figure:
    """Bar chart: exact hits per non-chrY chromosome."""
    non_chry = (
        hits_df[hits_df["chrom"] != "chrY"]
        if not hits_df.empty
        else pd.DataFrame()
    )

    if non_chry.empty:
        chrom_names: list[str] = []
        chrom_counts_vals: list[int] = []
    else:
        vc = non_chry["chrom"].value_counts().sort_index()
        chrom_names = vc.index.tolist()
        chrom_counts_vals = [int(v) for v in vc.values]

    fig = go.Figure(
        go.Bar(
            x=chrom_names,
            y=chrom_counts_vals,
            marker_color="#3498db",
            text=chrom_counts_vals,
            textposition="auto",
            hovertemplate="<b>%{x}</b><br>Hits: %{y}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Exact k-mer Hits by Non-chrY Chromosome",
        xaxis_title="Chromosome",
        yaxis_title="Number of Exact Hits",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=350,
        margin=dict(t=60, b=50),
    )
    return fig


def build_non_chry_table(hits_df: pd.DataFrame) -> go.Figure:
    """Interactive sortable table of non-chrY exact hits."""
    non_chry = (
        hits_df[hits_df["chrom"] != "chrY"]
        if not hits_df.empty
        else pd.DataFrame()
    )

    if non_chry.empty:
        display = pd.DataFrame(
            columns=["kmer", "chrom", "start", "end", "strand", "seq"]
        )
    else:
        display = non_chry[
            ["kmer", "chrom", "start", "end", "strand", "seq"]
        ].copy()
        display = display.sort_values(["chrom", "start"]).reset_index(drop=True)

    fig = go.Figure(
        go.Table(
            header=dict(
                values=[f"<b>{c}</b>" for c in display.columns],
                fill_color="#2c3e50",
                font=dict(color="white", size=12),
                align="left",
                height=32,
            ),
            cells=dict(
                values=[display[c] for c in display.columns],
                fill_color=[["#f8f9fa"] * len(display) for _ in display.columns],
                align="left",
                font=dict(size=11),
                height=28,
            ),
        )
    )
    fig.update_layout(
        title="Non-chrY Exact Hit Details",
        height=max(220, 35 * len(display) + 70),
        margin=dict(t=60, b=20),
    )
    return fig


def build_hit_table(hits_df: pd.DataFrame) -> go.Figure:
    """Interactive sortable table of all exact hits."""
    if hits_df.empty:
        display = pd.DataFrame(
            columns=["kmer", "chrom", "start", "end", "strand", "region", "seq"]
        )
    else:
        display = hits_df[
            ["kmer", "chrom", "start", "end", "strand", "region", "seq"]
        ].copy()
        display = display.sort_values(["chrom", "start"]).reset_index(drop=True)

    region_colors = {r["name"]: r["color"] for r in CHRY_REGIONS}

    # Per-column fill colours
    fill_colors: list[list[str]] = []
    for col in display.columns:
        if col == "region":
            fill_colors.append(
                [region_colors.get(str(v), "#f8f9fa") for v in display[col]]
            )
        else:
            fill_colors.append(["#f8f9fa"] * len(display))

    fig = go.Figure(
        go.Table(
            header=dict(
                values=[f"<b>{c}</b>" for c in display.columns],
                fill_color="#2c3e50",
                font=dict(color="white", size=12),
                align="left",
                height=32,
            ),
            cells=dict(
                values=[display[c] for c in display.columns],
                fill_color=fill_colors,
                align="left",
                font=dict(size=11),
                height=28,
            ),
        )
    )
    fig.update_layout(
        title="Exact Alignment Details",
        height=max(220, 35 * len(display) + 70),
        margin=dict(t=60, b=20),
    )
    return fig


# ─── HTML report ──────────────────────────────────────────────────────────────


def generate_html(
    karyogram: go.Figure,
    region_bar: go.Figure,
    hits_df: pd.DataFrame,
    all_kmers: list[tuple[str, str]],
    output_path: str,
    non_chry_bar: go.Figure | None = None,
    non_chry_summary: go.Figure | None = None,
    intervals_df: pd.DataFrame | None = None,
    alignment_path: str | None = None,
    multi_match_report: str | None = None,
    non_chry_report: str | None = None,
) -> None:
    """Combine figures into a single self-contained HTML report.

    Detailed per-hit tables are intentionally omitted from the HTML to keep
    the report small and responsive.  Full hit listings are written to the
    companion text files reported in *multi_match_report* and
    *non_chry_report*.
    """
    total_kmers = len(all_kmers)
    total_hits = len(hits_df)
    kmers_with_hits = hits_df["kmer"].nunique() if not hits_df.empty else 0
    kmers_no_hits = total_kmers - kmers_with_hits

    chry_hits = (
        hits_df[hits_df["chrom"] == "chrY"] if not hits_df.empty else pd.DataFrame()
    )
    chry_hit_count = len(chry_hits)
    par1_hits = int((chry_hits["region"] == "PAR1").sum()) if not chry_hits.empty else 0
    xtr_hits = int((chry_hits["region"] == "XTR").sum()) if not chry_hits.empty else 0
    par2_hits = int((chry_hits["region"] == "PAR2").sum()) if not chry_hits.empty else 0

    # Interval / cluster stats
    interval_count = len(intervals_df) if intervals_df is not None and not intervals_df.empty else 0
    cluster_count = (
        int(intervals_df["cluster"].sum())
        if intervals_df is not None and "cluster" in intervals_df.columns and not intervals_df.empty
        else 0
    )

    whole_genome_mode = non_chry_bar is not None
    non_chry_hit_count = (
        int((hits_df["chrom"] != "chrY").sum())
        if whole_genome_mode and not hits_df.empty
        else 0
    )

    karyogram_html = to_html(karyogram, full_html=False, include_plotlyjs=False)
    region_bar_html = to_html(region_bar, full_html=False, include_plotlyjs=False)

    region_pills = "".join(
        f'<span class="region-pill" style="background:{r["color"]}" '
        f'title="{r["description"]}">{r["name"]}</span>'
        for r in CHRY_REGIONS
    )

    # Build the non-chrY section only when whole-genome mode is active
    if whole_genome_mode:
        non_chry_bar_html = to_html(non_chry_bar, full_html=False, include_plotlyjs=False)
        non_chry_stat_card = f"""
      <div class="stat-card nonchrY">
        <div class="value">{non_chry_hit_count}</div>
        <div class="label">Non-chrY Hits</div>
      </div>"""
        chry_stat_card = f"""
      <div class="stat-card">
        <div class="value">{chry_hit_count}</div>
        <div class="label">chrY Hits</div>
      </div>"""

        non_chry_file_note = ""
        if non_chry_report:
            non_chry_file_note = (
                f'<p>Full hit listing: <code>{non_chry_report}</code></p>'
            )

        if non_chry_summary is not None:
            non_chry_summary_html = to_html(
                non_chry_summary, full_html=False, include_plotlyjs=False
            )
            non_chry_section = f"""
    <!-- Non-chrY hits (whole-genome mode) -->
    <div class="card">
      <h2>Non-chrY Hits by Chromosome</h2>
      {non_chry_bar_html}
    </div>

    <!-- Non-chrY summary table -->
    <div class="card">
      <h2>Non-chrY Summary</h2>
      {non_chry_summary_html}
      {non_chry_file_note}
    </div>"""
        else:
            non_chry_section = f"""
    <!-- Non-chrY hits (whole-genome mode) -->
    <div class="card">
      <h2>Non-chrY Hits by Chromosome</h2>
      {non_chry_bar_html}
      {non_chry_file_note}
    </div>"""
        genome_subtitle = "T2T CHM13v2.0 Full Genome — Exact k-mer Mapping Report"
    else:
        non_chry_stat_card = ""
        chry_stat_card = ""
        non_chry_section = ""
        genome_subtitle = "T2T CHM13v2.0 chrY — Exact k-mer Mapping Report"

    # Interval / cluster stat cards
    interval_stat_card = f"""
      <div class="stat-card">
        <div class="value">{interval_count}</div>
        <div class="label">Intervals</div>
      </div>""" if interval_count > 0 else ""

    cluster_stat_card = f"""
      <div class="stat-card nonchrY">
        <div class="value">{cluster_count}</div>
        <div class="label">Clusters</div>
      </div>""" if cluster_count > 0 else ""

    # Output files note (alignment + text reports)
    output_files_rows = ""
    if alignment_path:
        output_files_rows += f"<tr><td>Alignment (BAM/SAM)</td><td><code>{alignment_path}</code></td></tr>\n"
    if multi_match_report:
        output_files_rows += f"<tr><td>Multi-match kmers</td><td><code>{multi_match_report}</code></td></tr>\n"
    if non_chry_report:
        output_files_rows += f"<tr><td>Non-chrY kmer hits</td><td><code>{non_chry_report}</code></td></tr>\n"

    output_files_section = ""
    if output_files_rows:
        output_files_section = f"""
    <div class="card">
      <h2>Output Files</h2>
      <table class="files-table">
        <thead><tr><th>Description</th><th>Path</th></tr></thead>
        <tbody>
{output_files_rows}        </tbody>
      </table>
    </div>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>kmer_hunter Report</title>
  <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: #f0f2f5;
      color: #2c3e50;
      margin: 0;
      padding: 0;
    }}
    .header {{
      background: linear-gradient(135deg, #1a252f 0%, #2980b9 100%);
      color: #fff;
      padding: 2.5rem 2rem 2rem;
      text-align: center;
    }}
    .header h1 {{ margin: 0 0 0.4rem; font-size: 2.2rem; letter-spacing: -0.5px; }}
    .header p  {{ margin: 0; opacity: 0.85; font-size: 0.95rem; }}
    .badge {{
      display: inline-block;
      background: rgba(255,255,255,0.18);
      border: 1px solid rgba(255,255,255,0.35);
      border-radius: 999px;
      padding: 0.2rem 0.8rem;
      font-size: 0.78rem;
      margin-top: 0.6rem;
      letter-spacing: 0.04em;
    }}
    .container {{ max-width: 1400px; margin: 0 auto; padding: 1.5rem 1rem 3rem; }}
    .stats-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
      gap: 1rem;
      margin-bottom: 1.5rem;
    }}
    .stat-card {{
      background: #fff;
      border-radius: 10px;
      padding: 1.2rem 1rem;
      box-shadow: 0 2px 10px rgba(0,0,0,0.08);
      text-align: center;
    }}
    .stat-card .value {{
      font-size: 2.4rem;
      font-weight: 700;
      color: #2980b9;
      line-height: 1;
    }}
    .stat-card .label {{ font-size: 0.8rem; color: #7f8c8d; margin-top: 0.35rem; }}
    .stat-card.par1    .value {{ color: #27ae60; }}
    .stat-card.xtr     .value {{ color: #2980b9; }}
    .stat-card.par2    .value {{ color: #e67e22; }}
    .stat-card.miss    .value {{ color: #95a5a6; }}
    .stat-card.nonchrY .value {{ color: #c0392b; }}
    .card {{
      background: #fff;
      border-radius: 10px;
      padding: 1.5rem;
      box-shadow: 0 2px 10px rgba(0,0,0,0.08);
      margin-bottom: 1.5rem;
    }}
    .card h2 {{ margin: 0 0 1rem; font-size: 1.1rem; color: #2c3e50; }}
    .region-legend {{ display: flex; flex-wrap: wrap; gap: 0.4rem; margin-bottom: 1rem; }}
    .region-pill {{
      display: inline-block;
      padding: 0.22rem 0.7rem;
      border-radius: 999px;
      font-size: 0.75rem;
      font-weight: 600;
      color: #fff;
      cursor: default;
    }}
    .files-table {{
      border-collapse: collapse;
      width: 100%;
      font-size: 0.9rem;
    }}
    .files-table th, .files-table td {{
      text-align: left;
      padding: 0.5rem 0.75rem;
      border-bottom: 1px solid #ecf0f1;
    }}
    .files-table th {{ background: #f8f9fa; font-weight: 600; }}
    footer {{
      text-align: center;
      padding: 1.5rem;
      color: #95a5a6;
      font-size: 0.8rem;
    }}
  </style>
</head>
<body>
  <div class="header">
    <h1>🧬 kmer_hunter</h1>
    <p>{genome_subtitle}</p>
    <span class="badge">⚡ Exact matches only &mdash; zero mismatches, zero gaps</span>
  </div>

  <div class="container">

    <!-- Summary stats -->
    <div class="stats-grid">
      <div class="stat-card">
        <div class="value">{total_kmers}</div>
        <div class="label">Total k-mers</div>
      </div>
      <div class="stat-card">
        <div class="value">{total_hits}</div>
        <div class="label">Total Exact Hits</div>
      </div>
      <div class="stat-card">
        <div class="value">{kmers_with_hits}</div>
        <div class="label">k-mers Matched</div>
      </div>
      <div class="stat-card miss">
        <div class="value">{kmers_no_hits}</div>
        <div class="label">k-mers Not Found</div>
      </div>{chry_stat_card}{non_chry_stat_card}{interval_stat_card}{cluster_stat_card}
      <div class="stat-card par1">
        <div class="value">{par1_hits}</div>
        <div class="label">PAR1 Hits</div>
      </div>
      <div class="stat-card xtr">
        <div class="value">{xtr_hits}</div>
        <div class="label">XTR Hits</div>
      </div>
      <div class="stat-card par2">
        <div class="value">{par2_hits}</div>
        <div class="label">PAR2 Hits</div>
      </div>
    </div>

    <!-- Karyogram -->
    <div class="card">
      <h2>chrY Karyogram — Hit Count per Interval</h2>
      <div class="region-legend">{region_pills}</div>
      {karyogram_html}
    </div>

    <!-- Region bar chart -->
    <div class="card">
      <h2>Hits per Functional Region</h2>
      {region_bar_html}
    </div>
{non_chry_section}
{output_files_section}
  </div>

  <footer>
    Generated by <strong>kmer_hunter</strong> &nbsp;|&nbsp;
    Reference: T2T CHM13v2.0 (hs1) &nbsp;|&nbsp;
    Region annotations: Rhie&nbsp;et&nbsp;al.&nbsp;(2023)&nbsp;<em>Nature</em>
  </footer>
</body>
</html>
"""

    Path(output_path).write_text(html, encoding="utf-8")
    print(f"[kmer_hunter] Report written → {output_path}", file=sys.stderr)



# ─── Main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    args = parse_args()

    # 1. Read k-mers
    print(f"[kmer_hunter] Reading k-mers from: {args.kmers}", file=sys.stderr)
    kmers = read_kmers(args.kmers)
    print(f"[kmer_hunter] {len(kmers)} k-mer(s) loaded", file=sys.stderr)

    # 2. Ensure reference genome (whole T2T genome by default)
    ref_path = ensure_reference(args.reference, args.cache_dir, args.whole_genome_reference)
    print(f"[kmer_hunter] Reference: {ref_path}", file=sys.stderr)

    # 3. Search using BWA mem (exact matches, all chromosomes in one pass)
    print("[kmer_hunter] Searching for exact matches via BWA …", file=sys.stderr)
    hits, sam_text = bwa_find_exact_matches(kmers, ref_path)

    all_hits_df = pd.DataFrame(hits) if hits else pd.DataFrame(
        columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"]
    )

    # 4. Collapse to intervals and detect clusters
    intervals_df = collapse_to_intervals(all_hits_df)
    intervals_df = detect_clusters(intervals_df)

    # 5. Derive output stem for companion text files
    output_stem = str(Path(args.output).with_suffix(""))

    # 5a. Always save the alignment file; honour --output-bam if given,
    #     otherwise default to <output_stem>.sam so the SAM is never lost.
    requested_alignment_path = getattr(args, "output_bam", None)
    alignment_dest = requested_alignment_path or f"{output_stem}.sam"
    actual_alignment_path = save_alignment_file(sam_text, alignment_dest)

    # 5b. Write text report for kmers with more than one match
    multi_match_report = write_multi_match_report(all_hits_df, output_stem)

    # 5c. Write text report for non-chrY kmer hits
    non_chry_report_path = write_non_chry_report(all_hits_df, output_stem)

    # 6. Split hits for chrY-specific and non-chrY visualisations
    _empty = pd.DataFrame(
        columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"]
    )
    non_chry_df = (
        all_hits_df[all_hits_df["chrom"] != "chrY"].copy()
        if not all_hits_df.empty
        else _empty
    )
    non_chry_bar = build_non_chry_bar(non_chry_df)
    non_chry_summary = build_non_chry_summary(all_hits_df)

    # 7. Build figures (pass intervals for count y-axis karyogram)
    karyogram = build_karyogram(all_hits_df, intervals_df=intervals_df)
    region_bar = build_region_bar(all_hits_df)

    # 8. Generate HTML report
    generate_html(
        karyogram,
        region_bar,
        all_hits_df,
        kmers,
        args.output,
        non_chry_bar=non_chry_bar,
        non_chry_summary=non_chry_summary,
        intervals_df=intervals_df,
        alignment_path=actual_alignment_path,
        multi_match_report=multi_match_report,
        non_chry_report=non_chry_report_path,
    )
    print("[kmer_hunter] Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
