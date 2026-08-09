import gzip
import os
import tempfile
import unittest

from TEToolkit.GeneFeatures import GFF_Reader, GeneFeatures
from TEToolkit.TEindex import BinaryTree, TEfeatures


GENE_GTF = (
    'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "plus"\n'
    'chr1\ttest\texon\t150\t250\t.\t-\t.\tgene_id "minus"\n'
)
TE_GTF = (
    'chr1\ttest\texon\t1\t50000\t.\t+\t.\tgene_id "elem"; '
    'transcript_id "instance"; family_id "fam"; class_id "class"\n'
)


class AnnotationTests(unittest.TestCase):
    def test_plain_and_gzip_gtf_read_identically(self):
        with tempfile.TemporaryDirectory() as directory:
            plain = os.path.join(directory, "genes.gtf")
            compressed = plain + ".gz"
            with open(plain, "w") as handle:
                handle.write(GENE_GTF)
            with gzip.open(compressed, "wt") as handle:
                handle.write(GENE_GTF)
            self.assertEqual(list(GFF_Reader(plain, "gene_id")),
                             list(GFF_Reader(compressed, "gene_id")))

    def test_gene_annotation_honors_strand_and_inclusive_boundaries(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "genes.gtf")
            with open(filename, "w") as handle:
                handle.write(GENE_GTF)
            features = GeneFeatures(filename, "yes", "exon", "gene_id")
            self.assertEqual(features.Gene_annotation([["chr1", 100, 100, "+"]]), ['"plus"'])
            self.assertEqual(features.Gene_annotation([["chr1", 200, 200, "."]]),
                             ['"minus"', '"plus"'])

    def test_te_index_crosses_bins_groups_elements_and_finds_family(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "tes.gtf")
            with open(filename, "w") as handle:
                handle.write(TE_GTF)
            features = TEfeatures()
            features.build(filename)
            self.assertEqual(features.findOvpTE("chr1", 9999, 10001), [0])
            self.assertEqual(features.getFamilyID("chr1", 40000, 40001), "fam")
            self.assertEqual(features.groupByEle([2.5]), {"elem:fam:class": 2.5})

    def test_avl_rotations_keep_one_root_and_valid_parent_links(self):
        tree = BinaryTree()
        for index, start in enumerate((50001, 40001, 30001, 20001, 10001)):
            tree.insert(start, start + 10, index)
        nodes = []

        def visit(node):
            if node is None:
                return
            nodes.append(node)
            if node.left:
                self.assertIs(node.left.parent, node)
            if node.right:
                self.assertIs(node.right.parent, node)
            visit(node.left)
            visit(node.right)

        visit(tree.root)
        self.assertEqual(sum(node.isRoot() for node in nodes), 1)
        self.assertTrue(tree.root.isRoot())
        self.assertIn(tree.children_count(), (1, 2))


if __name__ == "__main__":
    unittest.main()
