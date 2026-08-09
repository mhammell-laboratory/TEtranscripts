import importlib.machinery
import importlib
import importlib.util
import os
import sys
import tempfile
import types
import unittest


REPOSITORY_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def load_script(name):
    path = os.path.join(REPOSITORY_ROOT, "bin", name)
    injected = False
    if "pysam" not in sys.modules:
        try:
            importlib.import_module("pysam")
        except ImportError:
            sys.modules["pysam"] = types.ModuleType("pysam")
            injected = True
    try:
        loader = importlib.machinery.SourceFileLoader("test_" + name, path)
        spec = importlib.util.spec_from_loader(loader.name, loader)
        module = importlib.util.module_from_spec(spec)
        loader.exec_module(module)
        return module
    finally:
        if injected:
            del sys.modules["pysam"]


class CountingCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modules = [load_script("TEcount"), load_script("TEtranscripts")]

    def test_cigar_coordinates_cover_every_reference_operation(self):
        cigar = [(4, 5), (0, 10), (1, 2), (7, 3), (2, 4),
                 (8, 2), (3, 20), (0, 5), (5, 3), (6, 2)]
        expected = [
            ["chr1", 100, 109, "+"], ["chr1", 110, 112, "+"],
            ["chr1", 117, 118, "+"], ["chr1", 139, 143, "+"],
        ]
        for module in self.modules:
            with self.subTest(module=module.__name__):
                self.assertEqual(module.fetch_exon("chr1", 99, cigar, 1, "BAM"), expected)

    def test_ambiguous_annotation_assignment_conserves_mass(self):
        for module in self.modules:
            with self.subTest(module=module.__name__):
                counts = {"A": 10.0, "B": 10.0}
                before = sum(counts.values())
                module.resolve_annotation_ambiguity(counts, [([["A"], ["B"]], 1.0)])
                self.assertAlmostEqual(sum(counts.values()) - before, 1.0)
                self.assertEqual(counts, {"A": 10.5, "B": 10.5})

                zero = {"A": 0.0, "B": 0.0}
                module.resolve_annotation_ambiguity(zero, [([["A", "B"]], 1.0)])
                self.assertEqual(zero, {"A": 0.5, "B": 0.5})

    def test_tecount_writes_a_rectangular_union_of_features(self):
        module = self.modules[0]
        with tempfile.TemporaryDirectory() as directory:
            prefix = os.path.join(directory, "counts")
            module.output_count_tbl({"s1": {"A": 1}, "s2": {"B": 2}}, prefix)
            with open(prefix + ".cntTable") as handle:
                rows = [line.rstrip().split("\t") for line in handle]
        self.assertEqual(rows[0], ["gene/TE", "s1", "s2"])
        self.assertEqual(rows[1:], [["A", "1", "0"], ["B", "0", "2"]])
        self.assertTrue(all(len(row) == 3 for row in rows))

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam is not installed")
    def test_tiny_bam_runs_through_annotation_and_evidence_weighted_em(self):
        import pysam
        from TEToolkit.GeneFeatures import GeneFeatures
        from TEToolkit.TEindex import TEfeatures

        with tempfile.TemporaryDirectory() as directory:
            gene_gtf = os.path.join(directory, "genes.gtf")
            te_gtf = os.path.join(directory, "tes.gtf")
            bam = os.path.join(directory, "reads.bam")
            with open(gene_gtf, "w") as handle:
                handle.write(
                    'chr1\ttest\texon\t700\t799\t.\t+\t.\tgene_id "gene1"\n'
                )
            with open(te_gtf, "w") as handle:
                handle.write(
                    'chr1\ttest\texon\t300\t399\t.\t+\t.\tgene_id "elem1"; '
                    'transcript_id "inst1"; family_id "fam"; class_id "class"\n'
                    'chr1\ttest\texon\t500\t599\t.\t+\t.\tgene_id "elem2"; '
                    'transcript_id "inst2"; family_id "fam"; class_id "class"\n'
                )

            header = {"HD": {"VN": "1.6", "SO": "queryname"},
                      "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(bam, "wb", header=header) as output:
                for name, start in (
                    ("a_unique_te", 299),
                    ("b_multi", 299),
                    ("b_multi", 499),
                    ("c_gene", 699),
                ):
                    read = pysam.AlignedSegment(output.header)
                    read.query_name = name
                    read.query_sequence = "A" * 50
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigartuples = [(0, 50)]
                    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                    output.write(read)

            for module in self.modules:
                with self.subTest(module=module.__name__):
                    genes = GeneFeatures(gene_gtf, "no", "exon", "gene_id")
                    tes = TEfeatures()
                    tes.build(te_gtf)
                    gene_counts, te_counts = module.count_transcript_abundance(
                        bam, "BAM", genes, tes, "no", "multi", False, 8, 50, 500
                    )
                    self.assertEqual(gene_counts['"gene1"'], 1.0)
                    self.assertAlmostEqual(sum(te_counts), 2.0)
                    self.assertGreater(te_counts[0], 1.5)
                    self.assertLess(te_counts[1], 0.5)


if __name__ == "__main__":
    unittest.main()
