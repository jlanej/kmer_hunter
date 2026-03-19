"""Unit tests for kmer_hunter.py.

All tests run offline — no network access or real reference downloads required.
"""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, call, patch

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
            hits = kh.bwa_find_exact_matches(kmers, str(ref))

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
            hits = kh.bwa_find_exact_matches([("k1", "ACGTACGTAC")], str(ref))

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

        with patch("subprocess.run", return_value=mock_fail):
            with self.assertRaises(SystemExit):
                kh.bwa_find_exact_matches([("k1", "ACGT")], str(ref))

    def test_empty_kmers_returns_empty(self):
        ref = self.tmpdir / "ref.fa"
        ref.write_text(">chrY\nACGT\n")
        Path(str(ref) + ".bwt").write_text("")

        mock_mem = MagicMock()
        mock_mem.returncode = 0
        mock_mem.stdout = "@HD\tVN:1.6\n"

        with patch("subprocess.run", return_value=mock_mem):
            hits = kh.bwa_find_exact_matches([], str(ref))

        self.assertEqual(hits, [])


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
                hits = kh.bwa_find_exact_matches([("kmer_1", target)], str(ref))

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
                hits = kh.bwa_find_exact_matches([("kmer_1", target)], str(ref))

        self.assertEqual(hits, [])


if __name__ == "__main__":
    unittest.main()
