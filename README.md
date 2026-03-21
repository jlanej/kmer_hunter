# kmer_hunter

Find exactly where k-mers live in the T2T CHM13v2.0 reference genome.

## Overview

`kmer_hunter` takes a list of k-mers, searches for **exact** (zero-mismatch,
zero-gap) occurrences in the T2T CHM13v2.0 chrY sequence, and generates a
visually rich, fully interactive HTML report containing:

- 🗺️ **chrY karyogram** — interactive Plotly chart with colour-coded PAR1,
  XTR, Ampliconic, Pericentromeric, Heterochromatin, PAR2 and Distal-Yq bands,
  with hit markers that show kmer name, position and region on hover.  Toggle
  between *Unique Hits Only* and *All Hits* views.
- 📊 **Region bar chart** — stacked bar chart of exact hit counts per
  functional region, with **Unique Hits** (k-mers mapping exactly once
  genome-wide, shown in the region colour) and **Multi-Hit** (k-mers with two
  or more genome-wide matches, shown in grey) distinguished by colour.
- 🔢 **Summary cards** — total k-mers, total hits, and unique hit counts for PAR1 / XTR / PAR2 (coloured to match the bar chart).

With `--whole-genome`, the report additionally includes:

- 🌐 **Non-chrY bar chart** — hit counts broken down by non-chrY chromosome,
  listed in natural genomic order (chr1, chr2, …, chr22, chrX, chrY, chrM).
- 📋 **Non-chrY sections** — a per-chromosome summary table plus a detailed
  non-chrY hit table with k-mer, chromosome, coordinates, strand and sequence.
- 🔢 **Non-chrY & chrY hit count cards** — at-a-glance comparison.

> **Exact matches only** — both forward and reverse-complement strands are
> searched; overlapping occurrences are all reported.  No mismatches, no gaps,
> no soft-clips.

---

## Quick start (Apptainer / Singularity)

```bash
apptainer run \
  --bind $PWD:/data \
  docker://ghcr.io/jlanej/kmer_hunter:latest \
  /data/kmers.txt \
  --output    /data/kmer_report.html \
  --cache-dir /data/.kmer_cache
```

Open `kmer_report.html` in any browser — no server required.

---

## Input format

**Plain text** — one k-mer per line:

```
ACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCA
GCTAGCTAGCTAGCTAGCTA
```

**FASTA** — standard FASTA:

```
>kmer_1
ACGTACGTACGTACGTACGT
>kmer_2
TGCATGCATGCATGCATGCA
```

---

## chrY region annotations

| Region | Colour | Coordinates (T2T hs1, approx.) | Description |
|---|---|---|---|
| PAR1 | 🟢 Green | 1 – 2,458,320 | Pseudoautosomal Region 1 (Yp telomere) |
| XTR | 🔵 Blue | 2,458,321 – 6,400,875 | X-Transposed Region |
| ↳ XTR1 | 🔵 Light blue | 2,727,073 – 5,914,561 | X-Transposed Region 1 |
| ↳ XTR2 | 🔵 Dark blue | 6,200,974 – 6,400,875 | X-Transposed Region 2 |
| Ampliconic | 🟣 Purple | 6,400,876 – 26,200,000 | Ampliconic sequences |
| Pericentromeric | 🔴 Red | 26,200,001 – 27,800,000 | Pericentromeric / Centromere |
| Heterochromatin | ⚫ Grey | 27,800,001 – 56,887,901 | Heterochromatin (DYZ1/DYZ2) |
| Distal Yq | 🔘 Light grey | 56,887,902 – 62,122,809 | Distal Yq (T2T-resolved) |
| PAR2 | 🟠 Orange | 62,122,810 – 62,460,029 | Pseudoautosomal Region 2 (Yq telomere) |

PAR coordinates from [chm13v2.0_PAR.bed](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed)
(T2T Consortium, Heng Li). The chrY lines of that BED file read:

```
chrY  0         2458320   (PAR1)
chrY  62122809  62460029  (PAR2)
```

XTR coordinates from Melissa Wilson (Arizona State Univ.) via
[GIAB genome-stratifications v3.1](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/CHM13v2.0/XY/T2T-CHM13v2.0_XY-stratifications.ipynb).
The released chrY XTR BED file is available at
[`CHM13v2.0_chrY_XTR.bed.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/CHM13v2.0/CHM13v2.0_chrY_XTR.bed.gz)
on the NCBI FTP.
The source file `CHM13XY_regions_from_Melissa.txt` defines two chrY XTR intervals:

```
HG002  chrY  2727072  5914561  XTR1
HG002  chrY  6200973  6400875  XTR2
```

K-mers falling inside one of these defined intervals are labelled **XTR1** or
**XTR2**; all other positions in the XTR span (2,458,321 – 6,400,875) are
labelled **XTR**.  In the summary statistics table, **XTR** shows the aggregate
count across the entire XTR region, with **XTR1** and **XTR2** broken out as
indented sub-rows.

Ampliconic, Pericentromeric, Heterochromatin and Distal Yq boundaries are
approximate, based on Rhie et al. (2023) *Nature* 621, 344–354, Fig. 1
([doi:10.1038/s41586-023-06457-y](https://doi.org/10.1038/s41586-023-06457-y)).

---

## Options

```
positional arguments:
  KMERS                 Input file: one k-mer per line or FASTA format

options:
  -o / --output PATH              Output HTML report [default: kmer_report.html]
  --reference PATH                Reference FASTA (.fa or .fa.gz). Defaults to T2T
                                  CHM13v2.0 chrY, auto-downloaded to --cache-dir.
  --cache-dir DIR                 Directory for cached reference genome
                                  [default: ~/.kmer_hunter_cache]
  --whole-genome                  Also search the full T2T CHM13v2.0 genome for
                                  non-chrY hits. Auto-downloads hs1.fa.gz (~830 MB)
                                  on first use unless --whole-genome-reference is set.
  --whole-genome-reference PATH   Path to a full-genome FASTA used with
                                  --whole-genome instead of auto-downloading.
```

---

## Whole-genome / non-chrY search

Pass `--whole-genome` to determine whether any of your k-mers also hit
chromosomes other than chrY.  This is useful for assessing chrY specificity.

```bash
apptainer run \
  --bind $PWD:/data \
  docker://ghcr.io/jlanej/kmer_hunter:latest \
  /data/kmers.txt \
  --output       /data/kmer_report.html \
  --cache-dir    /data/.kmer_cache \
  --whole-genome
```

On the first run this downloads the full T2T CHM13v2.0 genome (~830 MB
compressed) from UCSC to `--cache-dir`; subsequent runs use the cached file.
If you already have the genome locally, pass `--whole-genome-reference`:

```bash
kmer_hunter kmers.txt \
  --whole-genome \
  --whole-genome-reference /path/to/hs1.fa.gz \
  -o report.html
```

The report will include two additional sections:

1. **Non-chrY Hits by Chromosome** — bar chart of hit counts per chromosome
   (excluding chrY), listed in natural genomic order (chr1, chr2, …, chr22,
   chrX, chrY, chrM).
2. **Non-chrY Exact Hit Details** — full table of all hits outside chrY, also
   sorted in natural genomic chromosome order.

Summary stat cards for **chrY Hits** and **Non-chrY Hits** are also added so
the specificity of each k-mer is immediately visible.

---

## Reference genome caching

On first run, `kmer_hunter` automatically downloads the T2T CHM13v2.0 chrY
sequence (~16 MB compressed) from UCSC:

```
https://hgdownload.soe.ucsc.edu/goldenPath/hs1/chromosomes/chrY.fa.gz
```

When `--whole-genome` is used the full genome (~830 MB compressed) is cached:

```
https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
```

Subsequent runs re-use the cached file.  When running inside a container,
mount a host directory and pass `--cache-dir` so the download persists:

```bash
apptainer run \
  --bind $PWD:/data \
  docker://ghcr.io/jlanej/kmer_hunter:latest \
  /data/kmers.txt \
  --output    /data/kmer_report.html \
  --cache-dir /data/.kmer_cache
```

---

## Docker (without Apptainer)

```bash
docker run --rm \
  -v $PWD:/data \
  ghcr.io/jlanej/kmer_hunter:latest \
  /data/kmers.txt \
  --output    /data/kmer_report.html \
  --cache-dir /data/.kmer_cache
```

---

## Build from source

```bash
git clone https://github.com/jlanej/kmer_hunter
cd kmer_hunter
docker build -t kmer_hunter .
docker run --rm -v $PWD:/data kmer_hunter /data/kmers.txt -o /data/report.html
```

---

## References

- Rhie A. et al. (2023) "The complete sequence of a human Y chromosome"
  *Nature* 621, 344–354. <https://doi.org/10.1038/s41586-023-06457-y>
- T2T CHM13v2.0 (hs1): <https://github.com/marbl/CHM13>
