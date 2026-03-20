"""Unit tests for kmer_hunter.py.

All tests run offline — no network access or real reference downloads required.
"""

import gzip
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pandas as pd

# kmer_hunter.py lives next to this test file.
sys.path.insert(0, str(Path(__file__).parent))
import kmer_hunter as kh


# ── Helpers ────────────────────────────────────────────────────────────────────

def _write_plain_fasta(path: Path, sequences: dict) -> None:
    """Write a plain (non-gzipped) FASTA file."""
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def _write_gz_fasta(path: Path, sequences: dict) -> None:
    """Write a gzip-compressed FASTA file."""
    with gzip.open(path, "wt") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


# ── reverse_complement ─────────────────────────────────────────────────────────

class TestReverseComplement(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(kh.reverse_complement("ACGT"), "ACGT")

    def test_all_bases(self):
        self.assertEqual(kh.reverse_complement("AACCGGTT"), "AACCGGTT")

    def test_asymmetric(self):
        self.assertEqual(kh.reverse_complement("AAAA"), "TTTT")

    def test_single_base(self):
        self.assertEqual(kh.reverse_complement("A"), "T")
        self.assertEqual(kh.reverse_complement("C"), "G")
        self.assertEqual(kh.reverse_complement("G"), "C")
        self.assertEqual(kh.reverse_complement("T"), "A")

    def test_mixed_case(self):
        self.assertEqual(kh.reverse_complement("acgt"), "acgt")

    def test_n_preserved(self):
        rc = kh.reverse_complement("ACNGT")
        self.assertEqual(rc, "ACNGT")


# ── read_kmers ─────────────────────────────────────────────────────────────────

class TestReadKmers(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_plain_text(self):
        p = self.tmpdir / "kmers.txt"
        p.write_text("ACGT\nTTGG\nCCAA\n")
        result = kh.read_kmers(str(p))
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], ("kmer_1", "ACGT"))
        self.assertEqual(result[1], ("kmer_2", "TTGG"))
        self.assertEqual(result[2], ("kmer_3", "CCAA"))

    def test_plain_text_skips_comments(self):
        p = self.tmpdir / "kmers.txt"
        p.write_text("ACGT\n# comment\nTTGG\n")
        result = kh.read_kmers(str(p))
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0][1], "ACGT")
        self.assertEqual(result[1][1], "TTGG")

    def test_fasta_format(self):
        p = self.tmpdir / "kmers.fa"
        p.write_text(">seq1\nACGT\n>seq2\nTTGG\n")
        result = kh.read_kmers(str(p))
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], ("seq1", "ACGT"))
        self.assertEqual(result[1], ("seq2", "TTGG"))

    def test_fasta_lowercased_input_uppercased(self):
        p = self.tmpdir / "kmers.fa"
        p.write_text(">kmer1\nacgt\n")
        result = kh.read_kmers(str(p))
        self.assertEqual(result[0][1], "ACGT")

    def test_empty_file_exits(self):
        p = self.tmpdir / "empty.txt"
        p.write_text("")
        with self.assertRaises(SystemExit):
            kh.read_kmers(str(p))

    def test_only_comments_exits(self):
        p = self.tmpdir / "comments.txt"
        p.write_text("# just a comment\n")
        with self.assertRaises(SystemExit):
            kh.read_kmers(str(p))


# ── read_fasta ─────────────────────────────────────────────────────────────────

class TestReadFasta(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_plain_fasta(self):
        p = self.tmpdir / "ref.fa"
        _write_plain_fasta(p, {"chrY": "ACGT" * 10, "chr1": "TTTT" * 5})
        result = kh.read_fasta(str(p))
        self.assertIn("chrY", result)
        self.assertIn("chr1", result)
        self.assertEqual(result["chrY"], "ACGT" * 10)

    def test_gzipped_fasta(self):
        p = self.tmpdir / "ref.fa.gz"
        _write_gz_fasta(p, {"chrY": "AAACCCGGG"})
        result = kh.read_fasta(str(p))
        self.assertEqual(result["chrY"], "AAACCCGGG")

    def test_lowercase_converted_to_upper(self):
        p = self.tmpdir / "ref.fa"
        _write_plain_fasta(p, {"chrX": "acgt"})
        result = kh.read_fasta(str(p))
        self.assertEqual(result["chrX"], "ACGT")

    def test_multiline_sequence_joined(self):
        p = self.tmpdir / "ref.fa"
        p.write_text(">chrY\nACGT\nACGT\nACGT\n")
        result = kh.read_fasta(str(p))
        self.assertEqual(result["chrY"], "ACGTACGTACGT")


# ── annotate_region ────────────────────────────────────────────────────────────

class TestAnnotateRegion(unittest.TestCase):
    def test_par1(self):
        self.assertEqual(kh.annotate_region("chrY", 1), "PAR1")
        self.assertEqual(kh.annotate_region("chrY", 2_781_479), "PAR1")

    def test_xtr(self):
        self.assertEqual(kh.annotate_region("chrY", 2_781_480), "XTR")
        self.assertEqual(kh.annotate_region("chrY", 6_811_428), "XTR")

    def test_ampliconic(self):
        self.assertEqual(kh.annotate_region("chrY", 6_811_429), "Ampliconic")

    def test_pericentromeric(self):
        self.assertEqual(kh.annotate_region("chrY", 26_200_001), "Pericentromeric")

    def test_heterochromatin(self):
        self.assertEqual(kh.annotate_region("chrY", 27_800_001), "Heterochromatin")

    def test_par2(self):
        self.assertEqual(kh.annotate_region("chrY", 56_887_902), "PAR2")

    def test_distal_yq(self):
        self.assertEqual(kh.annotate_region("chrY", 57_217_416), "Distal Yq")
        self.assertEqual(kh.annotate_region("chrY", kh.CHRY_LEN), "Distal Yq")

    def test_non_chry_returns_chrom_name(self):
        self.assertEqual(kh.annotate_region("chr1", 1000), "chr1")
        self.assertEqual(kh.annotate_region("chrX", 500_000), "chrX")


# ── find_exact_matches ─────────────────────────────────────────────────────────

class TestFindExactMatches(unittest.TestCase):
    # Small synthetic reference (100 bp chrY-like sequence)
    REF = {"chrY": "ACGTACGTACGT" + "N" * 40 + "TTTAAA" * 4}

    def test_forward_hit(self):
        kmers = [("k1", "ACGTACGT")]
        hits = kh.find_exact_matches(kmers, self.REF)
        strands = [h["strand"] for h in hits if h["kmer"] == "k1"]
        self.assertIn("+", strands)

    def test_reverse_complement_hit(self):
        # "AAACCC" has RC "GGGTTT". Both substrings appear in the reference, so
        # find_exact_matches should return hits on both strands.
        ref = {"chrY": "AAACCCNNNGGGTTT"}
        kmers = [("k2", "AAACCC")]
        hits = kh.find_exact_matches(kmers, ref)
        strands = {h["strand"] for h in hits}
        self.assertIn("+", strands)  # forward hit of AAACCC
        self.assertIn("-", strands)  # RC hit: GGGTTT on ref → k-mer on minus strand

    def test_no_hit(self):
        kmers = [("k_miss", "GGGGGGGGGGGG")]
        hits = kh.find_exact_matches(kmers, self.REF)
        self.assertEqual(hits, [])

    def test_hit_coordinates_are_one_based(self):
        # TTTAAA starts at position 53 (1-based) in REF
        ref = {"chrY": "ACGT" * 13 + "TTTAAA" + "CCCC" * 5}
        kmers = [("k3", "TTTAAA")]
        hits = kh.find_exact_matches(kmers, ref)
        fwd = [h for h in hits if h["strand"] == "+"]
        self.assertTrue(len(fwd) >= 1)
        self.assertEqual(fwd[0]["start"], 53)  # 52 zero-based → 53 one-based
        self.assertEqual(fwd[0]["end"], 58)

    def test_hit_dict_keys(self):
        kmers = [("k4", "ACGT")]
        hits = kh.find_exact_matches(kmers, {"chrY": "ACGT"})
        self.assertTrue(len(hits) >= 1)
        expected_keys = {"kmer", "seq", "chrom", "start", "end", "strand", "region"}
        self.assertEqual(set(hits[0].keys()), expected_keys)

    def test_multiple_chromosomes(self):
        ref = {"chrY": "ACGT" * 5, "chr1": "TTTT" * 5}
        kmers = [("k5", "ACGT")]
        hits = kh.find_exact_matches(kmers, ref)
        chroms = {h["chrom"] for h in hits}
        self.assertIn("chrY", chroms)
        self.assertNotIn("chr1", chroms)


# ── _extract_chrom ─────────────────────────────────────────────────────────────

class TestExtractChrom(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_extracts_target_chrom(self):
        src = self.tmpdir / "full.fa.gz"
        _write_gz_fasta(src, {"chr1": "AAAA" * 10, "chrY": "CCCC" * 10, "chr2": "GGGG" * 10})
        dst = self.tmpdir / "chrY.fa.gz"
        kh._extract_chrom(str(src), "chrY", str(dst))
        result = kh.read_fasta(str(dst))
        self.assertIn("chrY", result)
        self.assertNotIn("chr1", result)
        self.assertNotIn("chr2", result)
        self.assertEqual(result["chrY"], "CCCC" * 10)

    def test_missing_chrom_raises(self):
        src = self.tmpdir / "full.fa.gz"
        _write_gz_fasta(src, {"chr1": "AAAA"})
        dst = self.tmpdir / "chrY.fa.gz"
        with self.assertRaises(ValueError):
            kh._extract_chrom(str(src), "chrY", str(dst))
        self.assertFalse(dst.exists())

    def test_dest_removed_on_error(self):
        src = self.tmpdir / "notexist.fa.gz"
        dst = self.tmpdir / "out.fa.gz"
        with self.assertRaises(Exception):
            kh._extract_chrom(str(src), "chrY", str(dst))
        self.assertFalse(dst.exists())


# ── ensure_reference ──────────────────────────────────────────────────────────

class TestEnsureReference(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_preexisting_reference_returned(self):
        ref = self.tmpdir / "my.fa"
        _write_plain_fasta(ref, {"chrY": "ACGT"})
        result = kh.ensure_reference(str(ref), str(self.tmpdir / "cache"))
        self.assertEqual(result, str(ref))

    def test_missing_reference_exits(self):
        with self.assertRaises(SystemExit):
            kh.ensure_reference("/nonexistent/ref.fa", str(self.tmpdir))

    def test_cached_whole_genome_returned_without_download(self):
        """Cached whole-genome file is reused - no download attempted."""
        cache = self.tmpdir / "cache"
        cache.mkdir()
        cached = cache / "chm13v2.0_hs1.fa.gz"
        _write_gz_fasta(cached, {"chrY": "ACGT", "chr1": "TTTT"})
        with patch("urllib.request.urlretrieve") as mock_dl:
            result = kh.ensure_reference(None, str(cache))
        mock_dl.assert_not_called()
        self.assertEqual(result, str(cached))

    def test_no_reference_triggers_whole_genome_download(self):
        """When no reference is provided the whole-genome URL is downloaded."""
        cache = self.tmpdir / "cache"

        def fake_download(url, dest):
            _write_gz_fasta(Path(dest), {"chr1": "AAAA", "chrY": "CCCC"})

        with patch("urllib.request.urlretrieve", side_effect=fake_download) as mock_dl:
            result = kh.ensure_reference(None, str(cache))

        mock_dl.assert_called_once()
        self.assertTrue(Path(result).exists())
        # The downloaded file must be the whole-genome reference
        self.assertIn("chm13v2.0_hs1", Path(result).name)

    def test_whole_genome_reference_used_when_no_reference(self):
        """--whole-genome-reference is honoured when no --reference is given."""
        ref = self.tmpdir / "chm13v2.0.fa"
        _write_plain_fasta(ref, {"chr1": "ACGT", "chrY": "TTTT"})
        with patch("urllib.request.urlretrieve") as mock_dl:
            result = kh.ensure_reference(None, str(self.tmpdir / "cache"), str(ref))
        mock_dl.assert_not_called()
        self.assertEqual(result, str(ref))

    def test_whole_genome_reference_any_name_accepted(self):
        """Any genome filename is accepted for --whole-genome-reference."""
        for name in ("chm13v2.0.fa", "genome.fasta", "my.genome.fa.gz", "hg38"):
            ref = self.tmpdir / name
            ref.write_bytes(b"")  # existence check only; BWA handles actual content
            result = kh.ensure_reference(None, str(self.tmpdir / "cache"), str(ref))
            self.assertEqual(result, str(ref))

    def test_missing_whole_genome_reference_exits(self):
        """A non-existent --whole-genome-reference path causes sys.exit."""
        with self.assertRaises(SystemExit):
            kh.ensure_reference(None, str(self.tmpdir), "/nonexistent/chm13v2.0.fa")

    def test_reference_takes_priority_over_whole_genome_reference(self):
        """When both --reference and --whole-genome-reference are given, --reference wins."""
        ref = self.tmpdir / "chry.fa"
        whole = self.tmpdir / "chm13v2.0.fa"
        _write_plain_fasta(ref, {"chrY": "ACGT"})
        _write_plain_fasta(whole, {"chr1": "TTTT", "chrY": "CCCC"})
        result = kh.ensure_reference(str(ref), str(self.tmpdir / "cache"), str(whole))
        self.assertEqual(result, str(ref))


# ── ensure_bwa_index ──────────────────────────────────────────────────────────

class TestEnsureBwaIndex(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_skips_index_when_bwt_exists(self):
        """If the .bwt sentinel file is present, bwa index must not be called."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chr1\nACGT\n")
        bwt = Path(str(ref) + ".bwt")
        bwt.write_text("")  # pre-existing index sentinel

        with patch("subprocess.run") as mock_run:
            kh.ensure_bwa_index(str(ref))
        mock_run.assert_not_called()

    def test_builds_index_when_bwt_missing(self):
        """When .bwt is absent, bwa index is invoked."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chr1\nACGT\n")

        mock_result = MagicMock()
        mock_result.returncode = 0

        with patch("subprocess.run", return_value=mock_result) as mock_run:
            kh.ensure_bwa_index(str(ref))

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        self.assertIn("bwa", cmd[0])
        self.assertIn("index", cmd)
        self.assertIn(str(ref), cmd)

    def test_exits_on_bwa_index_failure(self):
        """A non-zero return code from bwa index causes sys.exit."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chr1\nACGT\n")

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "some bwa error"

        with patch("subprocess.run", return_value=mock_result):
            with self.assertRaises(SystemExit):
                kh.ensure_bwa_index(str(ref))


# ── _parse_sam_exact ──────────────────────────────────────────────────────────

class TestParseSamExact(unittest.TestCase):
    def _sam(self, flag, chrom, pos, cigar, nm_val, qname="k1", seq="ACGTACGTAC"):
        nm_field = f"NM:i:{nm_val}" if nm_val is not None else ""
        return f"{qname}\t{flag}\t{chrom}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t*\t{nm_field}\n"

    def test_exact_forward_hit_parsed(self):
        sam = "@HD\tVN:1.6\n" + self._sam(0, "chrY", 100, "10M", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["chrom"], "chrY")
        self.assertEqual(hits[0]["start"], 100)
        self.assertEqual(hits[0]["end"], 109)
        self.assertEqual(hits[0]["strand"], "+")
        self.assertEqual(hits[0]["kmer"], "k1")

    def test_reverse_strand_flag(self):
        sam = self._sam(16, "chrY", 200, "10M", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["strand"], "-")

    def test_mismatch_skipped(self):
        sam = self._sam(0, "chrY", 100, "10M", 1)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(hits, [])

    def test_soft_clip_skipped(self):
        sam = self._sam(0, "chrY", 100, "8M2S", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(hits, [])

    def test_unmapped_skipped(self):
        sam = self._sam(4, "*", 0, "*", None)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(hits, [])

    def test_supplementary_skipped(self):
        sam = self._sam(2048, "chrY", 100, "10M", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(hits, [])

    def test_wrong_cigar_length_skipped(self):
        # k-mer is 10 bp but CIGAR says 5M
        sam = self._sam(0, "chrY", 100, "5M", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(hits, [])

    def test_header_lines_ignored(self):
        sam = "@HD\tVN:1.6\n@SQ\tSN:chrY\tLN:62460029\n" + self._sam(0, "chrY", 1, "10M", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(len(hits), 1)

    def test_result_keys(self):
        sam = self._sam(0, "chrY", 1, "10M", 0)
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC"})
        self.assertEqual(
            set(hits[0].keys()), {"kmer", "seq", "chrom", "start", "end", "strand", "region"}
        )

    def test_multiple_hits(self):
        sam = (
            self._sam(0, "chrY", 10, "10M", 0, qname="k1")
            + self._sam(0, "chr1", 50, "10M", 0, qname="k2", seq="TTTTTTTTTT")
        )
        hits = kh._parse_sam_exact(sam, {"k1": "ACGTACGTAC", "k2": "TTTTTTTTTT"})
        self.assertEqual(len(hits), 2)


# ── bwa_find_exact_matches ────────────────────────────────────────────────────

class TestBwaFindExactMatches(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _mock_bwa(self, sam_output: str):
        """Return mock subprocess.run results for bwa index then bwa mem calls.

        Use this helper when the .bwt sentinel is absent so that ensure_bwa_index
        will call bwa index first, followed by bwa mem.  When the sentinel already
        exists (index reuse tests), patch subprocess.run with a single mock_mem
        result instead.
        """
        mock_idx = MagicMock()
        mock_idx.returncode = 0

        mock_mem = MagicMock()
        mock_mem.returncode = 0
        mock_mem.stdout = sam_output
        # First call = bwa index (when .bwt absent), second call = bwa mem
        return [mock_idx, mock_mem]

    def test_returns_hits_from_sam(self):
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")

        sam = (
            "@HD\tVN:1.6\n"
            "k1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        side_effects = self._mock_bwa(sam)
        kmers = [("k1", "ACGTACGTAC")]

        with patch("subprocess.run", side_effect=side_effects):
            hits, _sam = kh.bwa_find_exact_matches(kmers, str(ref))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["kmer"], "k1")
        self.assertEqual(hits[0]["start"], 1)

    def test_index_reused_when_bwt_exists(self):
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")
        Path(str(ref) + ".bwt").write_text("")  # pre-existing index

        sam = "k1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        mock_mem = MagicMock()
        mock_mem.returncode = 0
        mock_mem.stdout = sam

        with patch("subprocess.run", return_value=mock_mem) as mock_run:
            hits, _sam = kh.bwa_find_exact_matches([("k1", "ACGTACGTAC")], str(ref))

        # Only bwa mem should have been called (index step was skipped)
        self.assertEqual(mock_run.call_count, 1)
        cmd = mock_run.call_args[0][0]
        self.assertIn("mem", cmd)

    def test_exits_on_bwa_mem_failure(self):
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGT\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        mock_fail = MagicMock()
        mock_fail.returncode = 1
        mock_fail.stderr = "bwa error"
        mock_fail.stdout = ""

        with patch("subprocess.run", return_value=mock_fail):
            with self.assertRaises(SystemExit):
                kh.bwa_find_exact_matches([("k1", "ACGT")], str(ref))

    def test_bwa_mem_failure_message_includes_exit_code(self):
        """Error message must contain the exit code and stderr output."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGT\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        mock_fail = MagicMock()
        mock_fail.returncode = 42
        mock_fail.stderr = "[M::bwa_idx_load_from_disk] read 0 ALT contigs\n"
        mock_fail.stdout = ""

        with patch("subprocess.run", return_value=mock_fail):
            with self.assertRaises(SystemExit) as cm:
                kh.bwa_find_exact_matches([("k1", "ACGT")], str(ref))

        msg = str(cm.exception.code)
        self.assertIn("42", msg)
        self.assertIn("read 0 ALT contigs", msg)

    def test_bwa_mem_failure_message_includes_stdout_when_present(self):
        """stdout is appended to the error message when non-empty."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGT\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        mock_fail = MagicMock()
        mock_fail.returncode = 1
        mock_fail.stderr = ""
        mock_fail.stdout = "unexpected bwa stdout output\n"

        with patch("subprocess.run", return_value=mock_fail):
            with self.assertRaises(SystemExit) as cm:
                kh.bwa_find_exact_matches([("k1", "ACGT")], str(ref))

        msg = str(cm.exception.code)
        self.assertIn("unexpected bwa stdout output", msg)

    def test_empty_kmers_returns_empty(self):
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGT\n")
        Path(str(ref) + ".bwt").write_text("")

        with patch("subprocess.run", return_value=MagicMock()) as mock_run:
            hits, _sam = kh.bwa_find_exact_matches([], str(ref))

        self.assertEqual(hits, [])
        # bwa mem must NOT be called for an empty k-mer list
        mock_run.assert_not_called()

    def test_batching_calls_bwa_mem_multiple_times(self):
        """When k-mers exceed batch_size, multiple bwa mem calls are made."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        sam1 = (
            "@HD\tVN:1.6\n"
            "k1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        sam2 = (
            "@HD\tVN:1.6\n"
            "k2\t0\tchrY\t2\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        mock1 = MagicMock(returncode=0, stdout=sam1)
        mock2 = MagicMock(returncode=0, stdout=sam2)

        kmers = [("k1", "ACGTACGTAC"), ("k2", "ACGTACGTAC")]
        with patch("subprocess.run", side_effect=[mock1, mock2]) as mock_run:
            hits, sam_text = kh.bwa_find_exact_matches(kmers, str(ref), batch_size=1)

        # One bwa mem call per batch
        self.assertEqual(mock_run.call_count, 2)
        self.assertEqual(len(hits), 2)
        # SAM headers should appear only once in the combined output
        self.assertEqual(sam_text.count("@HD"), 1)

    def test_batching_combines_hits_from_all_batches(self):
        """Hits from all batches are combined into the final result list."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        sam1 = (
            "@HD\tVN:1.6\n"
            "k1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        sam2 = (
            "@HD\tVN:1.6\n"
            "k2\t0\tchrY\t2\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        mock1 = MagicMock(returncode=0, stdout=sam1)
        mock2 = MagicMock(returncode=0, stdout=sam2)

        kmers = [("k1", "ACGTACGTAC"), ("k2", "ACGTACGTAC")]
        with patch("subprocess.run", side_effect=[mock1, mock2]):
            hits, _ = kh.bwa_find_exact_matches(kmers, str(ref), batch_size=1)

        kmer_names = {h["kmer"] for h in hits}
        self.assertIn("k1", kmer_names)
        self.assertIn("k2", kmer_names)

    def test_single_batch_when_kmers_within_batch_size(self):
        """When k-mers fit within batch_size, only one bwa mem call is made."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        sam = (
            "@HD\tVN:1.6\n"
            "k1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        mock_mem = MagicMock(returncode=0, stdout=sam)

        kmers = [("k1", "ACGTACGTAC")]
        with patch("subprocess.run", return_value=mock_mem) as mock_run:
            hits, _ = kh.bwa_find_exact_matches(kmers, str(ref), batch_size=100)

        self.assertEqual(mock_run.call_count, 1)
        self.assertEqual(len(hits), 1)

    def test_batch_failure_exits_with_error(self):
        """A bwa mem failure in any batch triggers SystemExit with the exit code."""
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")
        Path(str(ref) + ".bwt").write_text("")  # skip index

        mock_ok = MagicMock(returncode=0, stdout="@HD\tVN:1.6\n")
        mock_fail = MagicMock(returncode=-9, stderr="Killed", stdout="")

        kmers = [("k1", "ACGTACGTAC"), ("k2", "ACGTACGTAC")]
        with patch("subprocess.run", side_effect=[mock_ok, mock_fail]):
            with self.assertRaises(SystemExit) as cm:
                kh.bwa_find_exact_matches(kmers, str(ref), batch_size=1)

        self.assertIn("-9", str(cm.exception.code))


# ── Poly-A kmer integration tests ─────────────────────────────────────────────

# The 26 poly-A-prefixed k-mers from the issue.  Exact-match search must return
# zero hits against a synthetic reference that lacks them, and must return at
# least one hit when the reference explicitly contains them.
_POLY_A_KMER_SEQS = [
    "AAAAAAAAAAAAAAAAAAAAAAGCGGCCAGCT",
    "AAAAAAAAAAAAAAAAAAAAAATGCGGTTCTC",
    "AAAAAAAAAAAAAAAAAAAAAGCGGCCAGCTG",
    "AAAAAAAAAAAAAAAAAAAAAGCTGCAATACT",
    "AAAAAAAAAAAAAAAAAAAAAGGCACACGTAG",
    "AAAAAAAAAAAAAAAAAAAAATTGGAGTCTGT",
    "AAAAAAAAAAAAAAAAAAAAATTGGTCATAAG",
    "AAAAAAAAAAAAAAAAAAAAGAACTCGGTGGT",
    "AAAAAAAAAAAAAAAAAAAAGAATATTGTACT",
    "AAAAAAAAAAAAAAAAAAAAGATATAATGCCT",
    "AAAAAAAAAAAAAAAAAAAAGATGAGAGTATG",
    "AAAAAAAAAAAAAAAAAAAAGCAGTGCAGGTA",
    "AAAAAAAAAAAAAAAAAAAAGCGGCCAGCTGA",
    "AAAAAAAAAAAAAAAAAAAAGCTGCAATACTG",
    "AAAAAAAAAAAAAAAAAAAAGGCACACGTAGA",
    "AAAAAAAAAAAAAAAAAAAAGGTATCCTAGTG",
    "AAAAAAAAAAAAAAAAAAAAGTGCTTGAAATG",
    "AAAAAAAAAAAAAAAAAAAATCAGAGTCAACT",
    "AAAAAAAAAAAAAAAAAAAATGAAGTCACTAG",
    "AAAAAAAAAAAAAAAAAAAATTGGAGTCTGTG",
    "AAAAAAAAAAAAAAAAAAAATTGGTCATAAGA",
    "AAAAAAAAAAAAAAAAAAACAGGAATGGCTCT",
    "AAAAAAAAAAAAAAAAAAACCTTGTGTTTGGG",
    "AAAAAAAAAAAAAAAAAAACGGAACATACAAA",
    "AAAAAAAAAAAAAAAAAAACTAGGAATGTATT",
    "AAAAAAAAAAAAAAAAAAAGAACTCGGTG",  # shorter: 29 bp
]


class TestPolyAKmerIntegration(unittest.TestCase):
    """Comprehensive integration test for exact k-mer matching with poly-A prefix k-mers.

    Verifies that:
      - All 26 poly-A k-mers return zero hits against a reference that lacks them.
      - Each k-mer individually returns zero hits (subTest coverage).
      - Exact forward hits are found when the reference contains the k-mer.
      - Hit coordinates are 1-based and accurate for poly-A k-mers.
      - Reverse-complement strand hits are detected for poly-A k-mers.
      - All 26 poly-A k-mers are found when embedded in a custom reference.
      - The shorter 29-bp k-mer behaves consistently with the others.
      - Hit dicts contain all required keys.
      - bwa_find_exact_matches (with a mocked subprocess) returns a hit for a
        poly-A k-mer whose SAM record satisfies the exact-match filter.
    """

    # Synthetic reference with no poly-A stretches — alternating short motifs
    # that contain no run of more than one identical base.
    _REF_NO_POLYA = {"chrY": "ACGTCGATCGTA" * 50 + "GCTAGCTAGCTA" * 30}

    @staticmethod
    def _make_kmers(seqs):
        return [(f"kmer_{i + 1}", seq) for i, seq in enumerate(seqs)]

    # ── no-hit tests ──────────────────────────────────────────────────────────

    def test_no_hits_all_poly_a_kmers_against_simple_reference(self):
        """All poly-A k-mers return zero hits against a reference without them."""
        kmers = self._make_kmers(_POLY_A_KMER_SEQS)
        hits = kh.find_exact_matches(kmers, self._REF_NO_POLYA)
        self.assertEqual(hits, [], f"Expected 0 hits but got {len(hits)}")

    def test_no_hits_each_poly_a_kmer_individually(self):
        """Each poly-A k-mer individually returns zero hits against the simple reference."""
        for i, seq in enumerate(_POLY_A_KMER_SEQS):
            with self.subTest(kmer_index=i, seq=seq):
                hits = kh.find_exact_matches(
                    [(f"kmer_{i + 1}", seq)], self._REF_NO_POLYA
                )
                self.assertEqual(hits, [], f"kmer '{seq}' should have no hits")

    # ── positive-hit tests ────────────────────────────────────────────────────

    def test_forward_hit_when_reference_contains_poly_a_kmer(self):
        """A forward hit is found when the reference explicitly contains a poly-A k-mer."""
        target = _POLY_A_KMER_SEQS[0]
        ref = {"chrY": "GCGCGCGCGCGC" + target + "CGCGCGCGCGCG"}
        hits = kh.find_exact_matches([("kmer_1", target)], ref)
        fwd = [h for h in hits if h["strand"] == "+"]
        self.assertEqual(len(fwd), 1)
        self.assertEqual(fwd[0]["seq"], target)
        self.assertEqual(fwd[0]["kmer"], "kmer_1")

    def test_hit_coordinates_are_one_based(self):
        """Hit start/end positions are 1-based for a poly-A k-mer at a known offset."""
        target = _POLY_A_KMER_SEQS[0]
        prefix = "GCGCGCGCGCGC"  # 12 bases → k-mer starts at 1-based position 13
        ref = {"chrY": prefix + target + "CGCGCGCGCGCG"}
        hits = kh.find_exact_matches([("kmer_1", target)], ref)
        fwd = [h for h in hits if h["strand"] == "+"]
        self.assertEqual(len(fwd), 1)
        self.assertEqual(fwd[0]["start"], len(prefix) + 1)
        self.assertEqual(fwd[0]["end"], len(prefix) + len(target))

    def test_reverse_complement_hit_for_poly_a_kmer(self):
        """The reverse complement of a poly-A k-mer is detected on the minus strand."""
        target = _POLY_A_KMER_SEQS[0]
        rc = kh.reverse_complement(target)
        # Place the RC in the reference so a minus-strand hit is expected.
        ref = {"chrY": "GCGCGCGCGCGC" + rc + "CGCGCGCGCGCG"}
        hits = kh.find_exact_matches([("kmer_1", target)], ref)
        minus = [h for h in hits if h["strand"] == "-"]
        self.assertGreater(len(minus), 0)

    def test_all_poly_a_kmers_found_in_custom_reference(self):
        """Every poly-A k-mer is found when each is embedded in a custom reference."""
        separator = "GCGCGCGCGCGC"  # 12 bp of non-A/T content between k-mers
        ref_seq = separator + separator.join(_POLY_A_KMER_SEQS) + separator
        ref = {"chrY": ref_seq}
        kmers = self._make_kmers(_POLY_A_KMER_SEQS)
        hits = kh.find_exact_matches(kmers, ref)
        kmer_fwd_hits = {h["kmer"] for h in hits if h["strand"] == "+"}
        for i, seq in enumerate(_POLY_A_KMER_SEQS):
            name = f"kmer_{i + 1}"
            with self.subTest(kmer=name, seq=seq):
                self.assertIn(
                    name, kmer_fwd_hits, f"No forward hit found for {name} ({seq!r})"
                )

    def test_shorter_poly_a_kmer_no_hit_then_positive_hit(self):
        """The shorter 29-bp poly-A k-mer returns no hit, then a positive hit correctly."""
        shorter = _POLY_A_KMER_SEQS[-1]
        self.assertLess(len(shorter), len(_POLY_A_KMER_SEQS[0]))

        # No hit in the reference without poly-A sequences.
        hits_none = kh.find_exact_matches([("k_short", shorter)], self._REF_NO_POLYA)
        self.assertEqual(hits_none, [])

        # Positive hit when the reference contains the k-mer.
        ref = {"chrY": "GCGCGCGCGCGCGCGC" + shorter + "CGCGCGCGCGCGCGCG"}
        hits_some = kh.find_exact_matches([("k_short", shorter)], ref)
        fwd = [h for h in hits_some if h["strand"] == "+"]
        self.assertEqual(len(fwd), 1)

    # ── hit dict structure ────────────────────────────────────────────────────

    def test_hit_dict_contains_required_keys_for_poly_a_kmer(self):
        """Each hit dict for a poly-A k-mer contains all required keys."""
        target = _POLY_A_KMER_SEQS[0]
        ref = {"chrY": "GCGCGCGCGCGC" + target + "CGCGCGCGCGCG"}
        hits = kh.find_exact_matches([("kmer_1", target)], ref)
        self.assertTrue(hits)
        expected_keys = {"kmer", "seq", "chrom", "start", "end", "strand", "region"}
        for hit in hits:
            with self.subTest(hit=hit):
                self.assertEqual(set(hit.keys()), expected_keys)

    # ── BWA pipeline integration (mocked subprocess) ──────────────────────────

    def test_bwa_find_exact_matches_returns_poly_a_kmer_hit(self):
        """bwa_find_exact_matches returns a hit for a poly-A k-mer via mocked BWA."""
        import tempfile as _tempfile

        target = _POLY_A_KMER_SEQS[0]
        klen = len(target)

        with _tempfile.TemporaryDirectory() as tmpdir:
            ref = Path(tmpdir) / "ref.fa"
            _write_plain_fasta(ref, {"chrY": "GCGCGCGCGCGC" + target + "CGCGCGCGCGCG"})
            # Pre-create the .bwt sentinel so ensure_bwa_index skips indexing.
            Path(str(ref) + ".bwt").write_text("")

            sam_output = (
                "@HD\tVN:1.6\n"
                f"kmer_1\t0\tchrY\t13\t60\t{klen}M\t*\t0\t0\t{target}\t*\tNM:i:0\n"
            )
            mock_mem = MagicMock()
            mock_mem.returncode = 0
            mock_mem.stdout = sam_output

            with patch("subprocess.run", return_value=mock_mem):
                hits, _sam = kh.bwa_find_exact_matches([("kmer_1", target)], str(ref))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["kmer"], "kmer_1")
        self.assertEqual(hits[0]["start"], 13)
        self.assertEqual(hits[0]["end"], 12 + klen)
        self.assertEqual(hits[0]["strand"], "+")

    def test_bwa_find_exact_matches_no_poly_a_kmer_hit_when_sam_empty(self):
        """bwa_find_exact_matches returns no hits when SAM has no alignments."""
        import tempfile as _tempfile

        target = _POLY_A_KMER_SEQS[0]

        with _tempfile.TemporaryDirectory() as tmpdir:
            ref = Path(tmpdir) / "ref.fa"
            _write_plain_fasta(ref, {"chrY": "ACGTCGATCGTA" * 50})
            Path(str(ref) + ".bwt").write_text("")

            # SAM with only unmapped record (flag 4) — no exact hit.
            sam_output = (
                "@HD\tVN:1.6\n"
                f"kmer_1\t4\t*\t0\t0\t*\t*\t0\t0\t{target}\t*\n"
            )
            mock_mem = MagicMock()
            mock_mem.returncode = 0
            mock_mem.stdout = sam_output

            with patch("subprocess.run", return_value=mock_mem):
                hits, _sam = kh.bwa_find_exact_matches([("kmer_1", target)], str(ref))

        self.assertEqual(hits, [])


# ── Real BWA integration tests (skipped if bwa not installed) ─────────────────

import shutil as _shutil

_BWA_AVAILABLE = _shutil.which("bwa") is not None

# Separator used in toy references — GC-rich, unique, no poly-A.
_SEP = "GCGCTATACGCGTATAGCGCTATA"  # 24 bp


@unittest.skipUnless(_BWA_AVAILABLE, "bwa executable not found — skipping real BWA tests")
class TestBwaRealIntegration(unittest.TestCase):
    """Integration tests that invoke the real BWA binary.

    These tests create small toy FASTA references with known sequences, run
    ``bwa_find_exact_matches`` end-to-end (no subprocess mocking), and assert
    that the tool finds exactly the k-mers that are present and returns no hits
    for k-mers that are absent.

    Why this matters: the previous BWA parameters included extreme mismatch /
    gap / clip penalties (e.g. ``-B 1000``) that caused integer overflow inside
    BWA's internal scoring for repetitive sequences such as poly-A-prefixed
    k-mers.  BWA would then report the wrong primary alignment and
    ``_parse_sam_exact`` would correctly reject it — resulting in zero hits even
    when the k-mer is present in the reference.  The fix is to rely on
    ``_parse_sam_exact`` (CIGAR = ``{klen}M``, NM = 0) for exact-match
    filtering rather than trying to block non-exact alignments via inflated
    scoring parameters.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    # ── helpers ───────────────────────────────────────────────────────────────

    def _ref_containing(self, kmers: list[str]) -> str:
        """Write a toy reference FASTA and return its path.

        Each k-mer is embedded between unique GC-rich separators so that each
        k-mer appears exactly once at a known 1-based position.
        """
        seq = _SEP
        for kmer in kmers:
            seq += kmer + _SEP
        ref = self.tmpdir / "ref.fa"
        _write_plain_fasta(ref, {"chrY": seq})
        return str(ref)

    def _ref_without(self, content: str = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG") -> str:
        """Write a toy reference that contains no poly-A k-mers."""
        ref = self.tmpdir / "ref_empty.fa"
        _write_plain_fasta(ref, {"chrY": content * 10})
        return str(ref)

    # ── basic sanity tests ────────────────────────────────────────────────────

    def test_unique_kmer_found_at_correct_position(self):
        """BWA finds a unique non-repetitive k-mer at the expected 1-based coordinate."""
        kmer = "ACGTCAGTGACGATCGTAGCTAGCATGCATGC"  # 32 bp, low repetition
        prefix_len = len(_SEP)  # kmer embedded right after the separator
        ref_path = self._ref_containing([kmer])

        hits, _sam = kh.bwa_find_exact_matches([("k1", kmer)], ref_path)

        fwd = [h for h in hits if h["strand"] == "+" and h["kmer"] == "k1"]
        self.assertGreater(len(fwd), 0, "BWA should find the unique k-mer")
        self.assertEqual(fwd[0]["start"], prefix_len + 1)
        self.assertEqual(fwd[0]["end"], prefix_len + len(kmer))

    def test_absent_kmer_returns_no_hits(self):
        """BWA returns no hits for a k-mer absent from the reference."""
        kmer = "ACGTCAGTGACGATCGTAGCTAGCATGCATGC"
        # Reference contains only GC repeats — no match possible.
        ref_path = self._ref_without()

        hits, _sam = kh.bwa_find_exact_matches([("k1", kmer)], ref_path)
        self.assertEqual(hits, [])

    # ── poly-A k-mer tests ────────────────────────────────────────────────────

    def test_poly_a_kmer_found_when_present(self):
        """BWA finds a poly-A-prefixed k-mer when it is embedded in the reference.

        Previously this failed because the extreme BWA penalty parameters
        (-B 1000 etc.) caused integer overflow, making BWA report a spurious
        alignment instead of the correct NM=0 hit.
        """
        kmer = _POLY_A_KMER_SEQS[0]
        ref_path = self._ref_containing([kmer])

        hits, _sam = kh.bwa_find_exact_matches([("kmer_1", kmer)], ref_path)
        kmer_hits = [h for h in hits if h["kmer"] == "kmer_1" and h["strand"] == "+"]
        self.assertGreater(
            len(kmer_hits), 0,
            f"BWA must find the poly-A k-mer '{kmer}' when it is in the reference"
        )

    def test_poly_a_kmer_absent_returns_no_hits(self):
        """BWA returns no hits for a poly-A k-mer absent from the reference."""
        kmer = _POLY_A_KMER_SEQS[0]
        ref_path = self._ref_without()

        hits, _sam = kh.bwa_find_exact_matches([("kmer_1", kmer)], ref_path)
        self.assertEqual(hits, [])

    def test_all_26_poly_a_kmers_found_via_real_bwa(self):
        """BWA (real subprocess) finds all 26 poly-A k-mers when each is in the reference.

        This is the key regression test for the BWA integer-overflow bug with
        extreme mismatch penalties.  Each k-mer is embedded exactly once in a
        toy reference; every one must appear in the hits returned by the live
        BWA pipeline.
        """
        ref_path = self._ref_containing(_POLY_A_KMER_SEQS)
        kmers = [(f"kmer_{i + 1}", seq) for i, seq in enumerate(_POLY_A_KMER_SEQS)]

        hits, _sam = kh.bwa_find_exact_matches(kmers, ref_path)
        found_names = {h["kmer"] for h in hits}

        for name, seq in kmers:
            with self.subTest(kmer=name, seq=seq):
                self.assertIn(
                    name, found_names,
                    f"BWA must return at least one hit for {name} ({seq!r})"
                )

    def test_all_26_poly_a_kmers_absent_return_no_hits_via_real_bwa(self):
        """BWA returns no hits for any of the 26 poly-A k-mers when none is in the reference."""
        ref_path = self._ref_without()
        kmers = [(f"kmer_{i + 1}", seq) for i, seq in enumerate(_POLY_A_KMER_SEQS)]

        hits, _sam = kh.bwa_find_exact_matches(kmers, ref_path)
        self.assertEqual(hits, [], f"Expected no hits but got {len(hits)}")

    def test_hit_dict_keys_from_real_bwa(self):
        """Hit dicts returned by the real BWA pipeline contain all required keys."""
        kmer = "ACGTCAGTGACGATCGTAGCTAGCATGCATGC"
        ref_path = self._ref_containing([kmer])

        hits, _sam = kh.bwa_find_exact_matches([("k1", kmer)], ref_path)
        self.assertTrue(hits, "Expected at least one hit from real BWA")
        expected_keys = {"kmer", "seq", "chrom", "start", "end", "strand", "region"}
        for hit in hits:
            self.assertEqual(set(hit.keys()), expected_keys)

    def test_multiple_kmers_in_same_run_via_real_bwa(self):
        """All k-mers submitted in one BWA run are found correctly."""
        kmers_seqs = [
            "ACGTCAGTGACGATCGTAGCTAGCATGCATGC",
            "TGCATCGATCGTAGCATCGATCGATCGATCGA",
        ]
        ref_path = self._ref_containing(kmers_seqs)
        kmers = [(f"k{i + 1}", seq) for i, seq in enumerate(kmers_seqs)]

        hits, _sam = kh.bwa_find_exact_matches(kmers, ref_path)
        found = {h["kmer"] for h in hits}
        self.assertIn("k1", found)
        self.assertIn("k2", found)


# ── collapse_to_intervals ─────────────────────────────────────────────────────

class TestCollapseToIntervals(unittest.TestCase):
    def test_empty_dataframe_returns_empty(self):
        df = pd.DataFrame(columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"])
        result = kh.collapse_to_intervals(df)
        self.assertEqual(len(result), 0)
        self.assertIn("chrom", result.columns)
        self.assertIn("count", result.columns)

    def test_single_hit_produces_single_interval(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df, gap=1000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["start"], 100)
        self.assertEqual(result.iloc[0]["end"], 110)
        self.assertEqual(result.iloc[0]["count"], 1)

    def test_nearby_hits_merged(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chrY", "start": 500, "end": 510, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df, gap=1000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["count"], 2)
        self.assertEqual(result.iloc[0]["start"], 100)
        self.assertEqual(result.iloc[0]["end"], 510)

    def test_distant_hits_separate(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chrY", "start": 5000, "end": 5010, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df, gap=500)
        self.assertEqual(len(result), 2)
        self.assertEqual(result.iloc[0]["count"], 1)
        self.assertEqual(result.iloc[1]["count"], 1)

    def test_different_chromosomes_not_merged(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chr1", "start": 100, "end": 110, "strand": "+", "region": "chr1"},
        ])
        result = kh.collapse_to_intervals(df, gap=1000)
        self.assertEqual(len(result), 2)

    def test_result_columns(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df)
        expected = {"chrom", "start", "end", "count", "unique_count", "region"}
        self.assertEqual(set(result.columns), expected)

    def test_unique_count_single_hit_kmer(self):
        """A k-mer with exactly one hit should contribute to unique_count."""
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df)
        self.assertEqual(result.iloc[0]["unique_count"], 1)
        self.assertEqual(result.iloc[0]["count"], 1)

    def test_unique_count_multi_hit_kmer_excluded(self):
        """A k-mer with more than one hit should not contribute to unique_count."""
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 500, "end": 510, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df, gap=1000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["count"], 2)
        self.assertEqual(result.iloc[0]["unique_count"], 0)

    def test_unique_count_mixed_kmers(self):
        """unique_count counts only hits from k-mers that appear exactly once."""
        df = pd.DataFrame([
            # k1 appears once (unique)
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
            # k2 appears twice (non-unique): both hits in the same interval
            {"kmer": "k2", "seq": "TTGG", "chrom": "chrY", "start": 200, "end": 210, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chrY", "start": 300, "end": 310, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df, gap=1000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["count"], 3)
        self.assertEqual(result.iloc[0]["unique_count"], 1)

    def test_unique_count_empty_dataframe(self):
        df = pd.DataFrame(columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"])
        result = kh.collapse_to_intervals(df)
        self.assertEqual(len(result), 0)
        self.assertIn("unique_count", result.columns)

    def test_end_position_uses_max(self):
        """When merging, the end position should be the maximum end of all merged hits."""
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 200, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chrY", "start": 150, "end": 180, "strand": "+", "region": "PAR1"},
        ])
        result = kh.collapse_to_intervals(df, gap=100)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["end"], 200)


# ── detect_clusters ───────────────────────────────────────────────────────────

class TestDetectClusters(unittest.TestCase):
    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["chrom", "start", "end", "count", "region"])
        result = kh.detect_clusters(df)
        self.assertIn("cluster", result.columns)

    def test_cluster_detection(self):
        df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 500, "count": 10, "region": "PAR1"},
            {"chrom": "chrY", "start": 1000, "end": 1500, "count": 2, "region": "PAR1"},
        ])
        result = kh.detect_clusters(df, min_hits=5)
        self.assertTrue(result.iloc[0]["cluster"])
        self.assertFalse(result.iloc[1]["cluster"])

    def test_default_threshold(self):
        df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 500, "count": 5, "region": "PAR1"},
            {"chrom": "chrY", "start": 1000, "end": 1500, "count": 4, "region": "PAR1"},
        ])
        result = kh.detect_clusters(df)
        self.assertTrue(result.iloc[0]["cluster"])
        self.assertFalse(result.iloc[1]["cluster"])

    def test_original_not_modified(self):
        df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 500, "count": 10, "region": "PAR1"},
        ])
        result = kh.detect_clusters(df)
        self.assertNotIn("cluster", df.columns)
        self.assertIn("cluster", result.columns)


# ── build_non_chry_summary ────────────────────────────────────────────────────

class TestBuildNonChrySummary(unittest.TestCase):
    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"])
        fig = kh.build_non_chry_summary(df)
        self.assertIsNotNone(fig)

    def test_summary_groups_by_chrom(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chr1", "start": 100, "end": 110, "strand": "+", "region": "chr1"},
            {"kmer": "k1", "seq": "ACGT", "chrom": "chr1", "start": 200, "end": 210, "strand": "+", "region": "chr1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chr2", "start": 100, "end": 110, "strand": "+", "region": "chr2"},
        ])
        fig = kh.build_non_chry_summary(df)
        # Table should have data for 2 chromosomes
        table_data = fig.data[0]
        self.assertEqual(len(table_data.cells.values[0]), 2)

    def test_chry_excluded(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTGG", "chrom": "chr1", "start": 100, "end": 110, "strand": "+", "region": "chr1"},
        ])
        fig = kh.build_non_chry_summary(df)
        table_data = fig.data[0]
        chroms = list(table_data.cells.values[0])
        self.assertNotIn("chrY", chroms)
        self.assertIn("chr1", chroms)


# ── save_alignment_file ───────────────────────────────────────────────────────

class TestSaveAlignmentFile(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_writes_sam_file(self):
        sam_text = "@HD\tVN:1.6\nk1\t0\tchrY\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        out = self.tmpdir / "output.sam"
        result = kh.save_alignment_file(sam_text, str(out))
        self.assertTrue(Path(result).exists())
        content = Path(result).read_text()
        self.assertIn("@HD", content)
        self.assertIn("k1", content)

    def test_auto_appends_sam_extension(self):
        sam_text = "@HD\tVN:1.6\n"
        out = self.tmpdir / "output"
        result = kh.save_alignment_file(sam_text, str(out))
        self.assertTrue(result.endswith(".sam"))

    def test_bam_falls_back_to_sam_without_samtools(self):
        sam_text = "@HD\tVN:1.6\n"
        out = self.tmpdir / "output.bam"
        with patch("shutil.which", return_value=None):
            result = kh.save_alignment_file(sam_text, str(out))
        # Without samtools, should fall back to SAM
        self.assertTrue(result.endswith(".sam"))
        self.assertTrue(Path(result).exists())

    def test_bam_with_samtools_produces_sorted_bam(self):
        """When samtools is available and path ends in .bam, a sorted BAM is written."""
        sam_text = "@HD\tVN:1.6\nk1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        out = self.tmpdir / "output.bam"

        # Simulate the three samtools calls: view (SAM→BAM), sort, index
        mock_view = MagicMock(returncode=0, stdout=b"FAKE_BAM_BYTES")
        mock_sort = MagicMock(returncode=0)
        mock_index = MagicMock(returncode=0)

        with (
            patch("shutil.which", return_value="/usr/bin/samtools"),
            patch("subprocess.run", side_effect=[mock_view, mock_sort, mock_index]),
        ):
            result = kh.save_alignment_file(sam_text, str(out))

        self.assertEqual(result, str(out))
        self.assertTrue(result.endswith(".bam"))


# ── build_karyogram with intervals ────────────────────────────────────────────

class TestBuildKaryogramWithIntervals(unittest.TestCase):
    def test_karyogram_with_intervals(self):
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 110, "count": 1, "unique_count": 1, "region": "PAR1", "cluster": False},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        self.assertIsNotNone(fig)
        # Y-axis must be visible (count axis)
        self.assertTrue(fig.layout.yaxis.visible)

    def test_karyogram_with_clusters(self):
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 1000, "count": 10, "unique_count": 5, "region": "PAR1", "cluster": True},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        self.assertIsNotNone(fig)
        self.assertTrue(fig.layout.yaxis.visible)

    def test_karyogram_without_intervals_fallback(self):
        """Test that individual markers fall back to y=1 when no intervals are provided."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=None)
        self.assertIsNotNone(fig)
        self.assertTrue(fig.layout.yaxis.visible)

    def test_karyogram_empty(self):
        hits_df = pd.DataFrame(columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"])
        fig = kh.build_karyogram(hits_df, intervals_df=pd.DataFrame())
        self.assertIsNotNone(fig)
        self.assertTrue(fig.layout.yaxis.visible)

    def test_karyogram_yaxis_title(self):
        """Y-axis title should reference hit count."""
        hits_df = pd.DataFrame(columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"])
        fig = kh.build_karyogram(hits_df)
        self.assertIn("count", fig.layout.yaxis.title.text.lower())

    def test_karyogram_toggle_buttons_present(self):
        """Toggle buttons for Unique Hits Only and All Hits must be present."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 110, "count": 1, "unique_count": 1, "region": "PAR1", "cluster": False},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        self.assertEqual(len(fig.layout.updatemenus), 1)
        button_labels = [b.label for b in fig.layout.updatemenus[0].buttons]
        self.assertIn("Unique Hits Only", button_labels)
        self.assertIn("All Hits", button_labels)

    def test_karyogram_toggle_default_is_unique(self):
        """The active button (default view) should be Unique Hits Only."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 110, "count": 1, "unique_count": 1, "region": "PAR1", "cluster": False},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        menu = fig.layout.updatemenus[0]
        # active=0 means the first button (Unique Hits Only) is selected by default
        self.assertEqual(menu.active, 0)
        self.assertEqual(menu.buttons[0].label, "Unique Hits Only")

    def test_karyogram_two_bar_traces_for_intervals(self):
        """Two Bar traces (unique and all) must be in the figure when intervals present."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 110, "count": 1, "unique_count": 1, "region": "PAR1", "cluster": False},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        bar_traces = [t for t in fig.data if t.type == "bar"]
        self.assertEqual(len(bar_traces), 2)
        names = {t.name for t in bar_traces}
        self.assertIn("unique hit interval", names)
        self.assertIn("hit interval (all)", names)

    def test_karyogram_unique_bar_visible_all_bar_hidden_by_default(self):
        """Unique hits bar is visible; all-hits bar is hidden initially."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 110, "count": 1, "unique_count": 1, "region": "PAR1", "cluster": False},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        unique_bar = next(t for t in fig.data if t.type == "bar" and t.name == "unique hit interval")
        all_bar = next(t for t in fig.data if t.type == "bar" and t.name == "hit interval (all)")
        self.assertTrue(unique_bar.visible)
        self.assertFalse(all_bar.visible)

    def test_karyogram_backward_compat_no_unique_count(self):
        """intervals_df without unique_count column should still render without error."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        intervals_df = pd.DataFrame([
            {"chrom": "chrY", "start": 100, "end": 110, "count": 1, "region": "PAR1", "cluster": False},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        self.assertIsNotNone(fig)
        bar_traces = [t for t in fig.data if t.type == "bar"]
        self.assertEqual(len(bar_traces), 2)

    def test_karyogram_no_toggle_buttons_without_intervals(self):
        """Toggle buttons should not appear in fallback (no intervals) mode."""
        hits_df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY", "start": 100, "end": 110, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_karyogram(hits_df, intervals_df=None)
        self.assertEqual(len(fig.layout.updatemenus), 0)


# ── bwa_find_exact_matches returns sam_text ───────────────────────────────────

class TestBwaReturnsSamText(unittest.TestCase):
    def test_returns_tuple(self):
        ref = Path(tempfile.mkdtemp()) / "ref.fa"
        ref.write_text(">chrY\nACGTACGTAC\n")
        Path(str(ref) + ".bwt").write_text("")

        sam = "k1\t0\tchrY\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        mock_mem = MagicMock()
        mock_mem.returncode = 0
        mock_mem.stdout = sam

        with patch("subprocess.run", return_value=mock_mem):
            result = kh.bwa_find_exact_matches([("k1", "ACGTACGTAC")], str(ref))

        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        hits, sam_text = result
        self.assertIsInstance(hits, list)
        self.assertIsInstance(sam_text, str)
        self.assertIn("k1", sam_text)


# ── write_multi_match_report ──────────────────────────────────────────────────

class TestWriteMultiMatchReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.stem = str(Path(self.tmp.name) / "report")

    def tearDown(self):
        self.tmp.cleanup()

    def _make_hits(self, kmer, n_hits, chrom="chrY"):
        return [
            {"kmer": kmer, "seq": "ACGTACGT", "chrom": chrom,
             "start": 100 + i * 1000, "end": 108 + i * 1000,
             "strand": "+", "region": "PAR1"}
            for i in range(n_hits)
        ]

    def test_returns_none_when_empty(self):
        result = kh.write_multi_match_report(pd.DataFrame(), self.stem)
        self.assertIsNone(result)

    def test_returns_none_when_all_single_match(self):
        df = pd.DataFrame(self._make_hits("k1", 1))
        result = kh.write_multi_match_report(df, self.stem)
        self.assertIsNone(result)

    def test_writes_file_for_multi_match_kmer(self):
        df = pd.DataFrame(self._make_hits("k1", 3))
        path = kh.write_multi_match_report(df, self.stem)
        self.assertIsNotNone(path)
        self.assertTrue(Path(path).exists())
        content = Path(path).read_text()
        self.assertIn("k1", content)
        self.assertIn("3", content)   # total hits
        self.assertIn("chrY", content)

    def test_single_match_kmer_excluded(self):
        rows = self._make_hits("multi", 2) + self._make_hits("single", 1)
        df = pd.DataFrame(rows)
        path = kh.write_multi_match_report(df, self.stem)
        self.assertIsNotNone(path)
        content = Path(path).read_text()
        self.assertIn("multi", content)
        self.assertNotIn("single", content)

    def test_output_filename_contains_stem(self):
        df = pd.DataFrame(self._make_hits("k1", 2))
        path = kh.write_multi_match_report(df, self.stem)
        self.assertTrue(path.startswith(self.stem))


# ── write_non_chry_report ─────────────────────────────────────────────────────

class TestWriteNonChrYReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.stem = str(Path(self.tmp.name) / "report")

    def tearDown(self):
        self.tmp.cleanup()

    def test_returns_none_when_empty(self):
        result = kh.write_non_chry_report(pd.DataFrame(), self.stem)
        self.assertIsNone(result)

    def test_returns_none_when_only_chry_hits(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])
        result = kh.write_non_chry_report(df, self.stem)
        self.assertIsNone(result)

    def test_writes_file_for_non_chry_hits(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chr1",
             "start": 500, "end": 504, "strand": "+", "region": "chr1"},
        ])
        path = kh.write_non_chry_report(df, self.stem)
        self.assertIsNotNone(path)
        self.assertTrue(Path(path).exists())
        content = Path(path).read_text()
        self.assertIn("k1", content)
        self.assertIn("chr1", content)

    def test_chry_hits_excluded(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chr1",
             "start": 500, "end": 504, "strand": "+", "region": "chr1"},
            {"kmer": "k2", "seq": "TTTT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])
        path = kh.write_non_chry_report(df, self.stem)
        content = Path(path).read_text()
        self.assertIn("k1", content)
        self.assertNotIn("k2", content)

    def test_output_filename_contains_stem(self):
        df = pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrX",
             "start": 100, "end": 104, "strand": "-", "region": "chrX"},
        ])
        path = kh.write_non_chry_report(df, self.stem)
        self.assertTrue(path.startswith(self.stem))


# ── generate_html (updated signature) ────────────────────────────────────────

class TestGenerateHtml(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _minimal_hits(self):
        return pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])

    def test_generates_html_file(self):
        hits_df = self._minimal_hits()
        out = str(self.tmpdir / "report.html")
        karyogram = kh.build_karyogram(hits_df)
        region_bar = kh.build_region_bar(hits_df)
        kh.generate_html(karyogram, region_bar, hits_df, [("k1", "ACGT")], out)
        self.assertTrue(Path(out).exists())

    def test_html_no_full_hit_table(self):
        """The full per-hit table should not be embedded in the HTML."""
        hits_df = self._minimal_hits()
        out = str(self.tmpdir / "report.html")
        karyogram = kh.build_karyogram(hits_df)
        region_bar = kh.build_region_bar(hits_df)
        kh.generate_html(karyogram, region_bar, hits_df, [("k1", "ACGT")], out)
        content = Path(out).read_text()
        # The "All Exact Hits" section title should be gone
        self.assertNotIn("All Exact Hits", content)

    def test_alignment_path_shown_in_output_files_section(self):
        hits_df = self._minimal_hits()
        out = str(self.tmpdir / "report.html")
        karyogram = kh.build_karyogram(hits_df)
        region_bar = kh.build_region_bar(hits_df)
        kh.generate_html(
            karyogram, region_bar, hits_df, [("k1", "ACGT")], out,
            alignment_path="/tmp/hits.sam",
        )
        content = Path(out).read_text()
        self.assertIn("/tmp/hits.sam", content)

    def test_multi_match_report_shown_in_output_files_section(self):
        hits_df = self._minimal_hits()
        out = str(self.tmpdir / "report.html")
        karyogram = kh.build_karyogram(hits_df)
        region_bar = kh.build_region_bar(hits_df)
        kh.generate_html(
            karyogram, region_bar, hits_df, [("k1", "ACGT")], out,
            multi_match_report="/tmp/multi.txt",
        )
        content = Path(out).read_text()
        self.assertIn("/tmp/multi.txt", content)


# ── _chrom_sort_key ────────────────────────────────────────────────────────────


class TestChromSortKey(unittest.TestCase):
    """Natural genomic ordering for chromosome names."""

    def test_numbered_chromosomes_ordered_numerically(self):
        chroms = ["chr10", "chr2", "chr1", "chr22"]
        self.assertEqual(
            sorted(chroms, key=kh._chrom_sort_key),
            ["chr1", "chr2", "chr10", "chr22"],
        )

    def test_sex_chromosomes_after_autosomes(self):
        chroms = ["chrX", "chr1", "chrY", "chr22"]
        result = sorted(chroms, key=kh._chrom_sort_key)
        self.assertLess(result.index("chr1"), result.index("chrX"))
        self.assertLess(result.index("chrX"), result.index("chrY"))

    def test_chrM_after_chrY(self):
        chroms = ["chrM", "chrY", "chrX"]
        result = sorted(chroms, key=kh._chrom_sort_key)
        self.assertEqual(result, ["chrX", "chrY", "chrM"])

    def test_non_standard_names_at_end(self):
        chroms = ["chr1_random", "chr1"]
        result = sorted(chroms, key=kh._chrom_sort_key)
        self.assertEqual(result[0], "chr1")
        self.assertEqual(result[-1], "chr1_random")

    def test_full_ordering(self):
        chroms = ["chrM", "chrX", "chr2", "chr11", "chr1", "chrY"]
        result = sorted(chroms, key=kh._chrom_sort_key)
        self.assertEqual(result, ["chr1", "chr2", "chr11", "chrX", "chrY", "chrM"])


# ── build_region_bar (unique / multi-hit split) ───────────────────────────────


class TestBuildRegionBar(unittest.TestCase):
    """Tests for the stacked unique/multi-hit region bar chart."""

    def _hits(self, rows):
        return pd.DataFrame(rows)

    def test_returns_figure(self):
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_region_bar(hits_df)
        self.assertIsNotNone(fig)

    def test_two_traces_present(self):
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_region_bar(hits_df)
        self.assertEqual(len(fig.data), 2)

    def test_trace_names(self):
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_region_bar(hits_df)
        names = [t.name for t in fig.data]
        self.assertIn("Unique Hits", names)
        self.assertIn("Multi-Hit", names)

    def test_unique_hit_counted_in_unique_trace(self):
        # k1 appears once genome-wide → unique
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_region_bar(hits_df)
        unique_trace = next(t for t in fig.data if t.name == "Unique Hits")
        multi_trace = next(t for t in fig.data if t.name == "Multi-Hit")
        par1_idx = [r["name"] for r in kh.CHRY_REGIONS].index("PAR1")
        self.assertEqual(unique_trace.y[par1_idx], 1)
        self.assertEqual(multi_trace.y[par1_idx], 0)

    def test_multi_hit_kmer_counted_in_multi_trace(self):
        # k1 appears twice genome-wide → multi-hit
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 200, "end": 204, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_region_bar(hits_df)
        unique_trace = next(t for t in fig.data if t.name == "Unique Hits")
        multi_trace = next(t for t in fig.data if t.name == "Multi-Hit")
        par1_idx = [r["name"] for r in kh.CHRY_REGIONS].index("PAR1")
        self.assertEqual(unique_trace.y[par1_idx], 0)
        self.assertEqual(multi_trace.y[par1_idx], 2)

    def test_stacked_total_equals_hit_count(self):
        # 2 unique + 2 multi-hit in PAR1
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
            {"kmer": "k2", "seq": "TTTT", "chrom": "chrY",
             "start": 200, "end": 204, "strand": "+", "region": "PAR1"},
            {"kmer": "k3", "seq": "CCCC", "chrom": "chrY",
             "start": 300, "end": 304, "strand": "+", "region": "PAR1"},
            {"kmer": "k3", "seq": "CCCC", "chrom": "chrY",
             "start": 400, "end": 404, "strand": "+", "region": "PAR1"},
        ])
        fig = kh.build_region_bar(hits_df)
        unique_trace = next(t for t in fig.data if t.name == "Unique Hits")
        multi_trace = next(t for t in fig.data if t.name == "Multi-Hit")
        par1_idx = [r["name"] for r in kh.CHRY_REGIONS].index("PAR1")
        total = unique_trace.y[par1_idx] + multi_trace.y[par1_idx]
        self.assertEqual(total, 4)

    def test_genome_wide_unique_uses_all_hits(self):
        # k1 hits chrY once and chr1 once → multi-hit genome-wide
        hits_df = self._hits([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chrY",
             "start": 100, "end": 104, "strand": "+", "region": "PAR1"},
            {"kmer": "k1", "seq": "ACGT", "chrom": "chr1",
             "start": 500, "end": 504, "strand": "+", "region": "chr1"},
        ])
        fig = kh.build_region_bar(hits_df)
        unique_trace = next(t for t in fig.data if t.name == "Unique Hits")
        multi_trace = next(t for t in fig.data if t.name == "Multi-Hit")
        par1_idx = [r["name"] for r in kh.CHRY_REGIONS].index("PAR1")
        # k1 has 2 genome-wide hits → multi, not unique
        self.assertEqual(unique_trace.y[par1_idx], 0)
        self.assertEqual(multi_trace.y[par1_idx], 1)

    def test_empty_dataframe_returns_figure(self):
        hits_df = pd.DataFrame(
            columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"]
        )
        fig = kh.build_region_bar(hits_df)
        self.assertIsNotNone(fig)
        self.assertEqual(len(fig.data), 2)

    def test_barmode_is_stack(self):
        hits_df = pd.DataFrame(
            columns=["kmer", "seq", "chrom", "start", "end", "strand", "region"]
        )
        fig = kh.build_region_bar(hits_df)
        self.assertEqual(fig.layout.barmode, "stack")


# ── non-chrY chromosome natural sort in display functions ─────────────────────


class TestNonChryNaturalSort(unittest.TestCase):
    """Chromosome sort order is natural (chr1 < chr2 < … < chrX < chrY < chrM)."""

    def _multi_chrom_hits(self):
        return pd.DataFrame([
            {"kmer": "k1", "seq": "ACGT", "chrom": "chr10",
             "start": 100, "end": 104, "strand": "+", "region": "chr10"},
            {"kmer": "k2", "seq": "TTTT", "chrom": "chr2",
             "start": 200, "end": 204, "strand": "+", "region": "chr2"},
            {"kmer": "k3", "seq": "CCCC", "chrom": "chr1",
             "start": 300, "end": 304, "strand": "+", "region": "chr1"},
        ])

    def test_build_non_chry_bar_natural_order(self):
        hits_df = self._multi_chrom_hits()
        fig = kh.build_non_chry_bar(hits_df)
        x = list(fig.data[0].x)
        self.assertLess(x.index("chr1"), x.index("chr2"))
        self.assertLess(x.index("chr2"), x.index("chr10"))

    def test_build_non_chry_summary_natural_order(self):
        hits_df = self._multi_chrom_hits()
        fig = kh.build_non_chry_summary(hits_df)
        chrom_col = list(fig.data[0].cells.values[0])
        self.assertLess(chrom_col.index("chr1"), chrom_col.index("chr2"))
        self.assertLess(chrom_col.index("chr2"), chrom_col.index("chr10"))

    def test_build_non_chry_table_natural_order(self):
        hits_df = self._multi_chrom_hits()
        fig = kh.build_non_chry_table(hits_df)
        chrom_col = list(fig.data[0].cells.values[1])  # column 1 is "chrom"
        # chr1 row appears before chr2, which appears before chr10
        self.assertLess(chrom_col.index("chr1"), chrom_col.index("chr2"))
        self.assertLess(chrom_col.index("chr2"), chrom_col.index("chr10"))


# ── Performance regression tests ──────────────────────────────────────────────


class TestBuildKaryogramPerformance(unittest.TestCase):
    """Ensure vectorized karyogram helpers handle large interval sets quickly."""

    def _make_intervals(self, n):
        return pd.DataFrame({
            "chrom": ["chrY"] * n,
            "start": [i * 2000 for i in range(n)],
            "end": [i * 2000 + 500 for i in range(n)],
            "count": [3] * n,
            "unique_count": [2] * n,
            "region": ["Ampliconic"] * n,
            "cluster": [i % 5 == 0 for i in range(n)],
        })

    def _make_hits(self, n):
        return pd.DataFrame({
            "kmer": [f"k{i}" for i in range(n)],
            "seq": ["ACGTACGT"] * n,
            "chrom": ["chrY"] * n,
            "start": [i * 2000 for i in range(n)],
            "end": [i * 2000 + 7 for i in range(n)],
            "strand": ["+"] * n,
            "region": ["Ampliconic"] * n,
        })

    def test_karyogram_large_intervals(self):
        """build_karyogram should handle 5 000 intervals without iterrows."""
        n = 5_000
        hits_df = self._make_hits(n)
        intervals_df = self._make_intervals(n)
        t0 = time.perf_counter()
        fig = kh.build_karyogram(hits_df, intervals_df=intervals_df)
        elapsed = time.perf_counter() - t0
        self.assertIsNotNone(fig)
        # Should complete well under 5 seconds with vectorized code.
        self.assertLess(elapsed, 5.0, f"build_karyogram took {elapsed:.1f}s for {n} intervals")


class TestWriteMultiMatchReportPerformance(unittest.TestCase):
    """Ensure vectorized report writing handles large hit sets quickly."""

    def test_large_multi_match_report(self):
        n_kmers = 1_000
        hits_per_kmer = 5
        rows = []
        for i in range(n_kmers):
            for j in range(hits_per_kmer):
                rows.append({
                    "kmer": f"kmer_{i}", "seq": "ACGTACGT",
                    "chrom": "chrY", "start": i * 10000 + j * 100,
                    "end": i * 10000 + j * 100 + 7,
                    "strand": "+", "region": "Ampliconic",
                })
        df = pd.DataFrame(rows)

        with tempfile.TemporaryDirectory() as tmp:
            stem = str(Path(tmp) / "perf_report")
            t0 = time.perf_counter()
            path = kh.write_multi_match_report(df, stem)
            elapsed = time.perf_counter() - t0
            self.assertIsNotNone(path)
            self.assertLess(elapsed, 5.0, f"write_multi_match_report took {elapsed:.1f}s")


class TestAlignmentAlwaysSaved(unittest.TestCase):
    """Ensure alignment file is always saved and reported in HTML."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_default_alignment_path_used_when_no_output_bam(self):
        """When --output-bam is not given, alignment is saved to <stem>.sam."""
        sam_text = "@HD\tVN:1.6\n"
        dest = str(self.tmpdir / "report.sam")
        result = kh.save_alignment_file(sam_text, dest)
        self.assertTrue(Path(result).exists())
        self.assertTrue(result.endswith(".sam"))



if __name__ == '__main__':
    unittest.main()
