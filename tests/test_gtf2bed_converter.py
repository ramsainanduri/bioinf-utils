import unittest
import os
import sys


# Add the scripts directory to the Python path
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts"))
)

from gtf2bed_converter import gtf_to_bed


class TestGTFToBED(unittest.TestCase):
    def setUp(self):
        # Create a sample GTF file for testing
        self.sample_gtf = "test.gtf"
        self.output_bed = "test_output.bed"
        with open(self.sample_gtf, "w") as f:
            f.write(
                'chr1\tsource\tgene\t100\t200\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1"; gene_type "protein_coding";\n'
                'chr1\tsource\tgene\t300\t400\t.\t+\t.\tgene_id "gene2"; gene_name "GENE2"; gene_type "lincRNA";\n'
                'chr1\tsource\texon\t100\t150\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1"; exon_number "1";\n'
                'chr1\tsource\texon\t150\t200\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1"; exon_number "2";\n'
                'chr1\tsource\tgene\t500\t600\t.\t+\t.\tgene_id "gene3"; gene_name "GENE3"; gene_type "pseudogene";\n'
            )

    def tearDown(self):
        # Remove the test files
        if os.path.exists(self.sample_gtf):
            os.remove(self.sample_gtf)
        if os.path.exists(self.output_bed):
            os.remove(self.output_bed)

    def test_simple_bed(self):
        """Test simple BED format conversion for genes."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            bed_columns="simple",
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0].strip(), "chr1\t99\t200\tGENE1")
        self.assertEqual(lines[1].strip(), "chr1\t299\t400\tGENE2")
        self.assertEqual(lines[2].strip(), "chr1\t499\t600\tGENE3")

    def test_detailed_bed(self):
        """Test detailed BED format conversion for genes."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            bed_columns="detailed",
            filter_key="gene_type",
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0].strip(), "chr1\t99\t200\t+\tGENE1\tprotein_coding")
        self.assertEqual(lines[1].strip(), "chr1\t299\t400\t+\tGENE2\tlincRNA")
        self.assertEqual(lines[2].strip(), "chr1\t499\t600\t+\tGENE3\tpseudogene")

    def test_include_patterns(self):
        """Test filtering with include patterns."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            filter_key="gene_type",
            include_patterns=["protein_coding"],
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].strip(), "chr1\t99\t200\tGENE1")

    def test_exclude_patterns(self):
        """Test filtering with exclude patterns."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            filter_key="gene_type",
            exclude_patterns=["protein_coding"],
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0].strip(), "chr1\t299\t400\tGENE2")
        self.assertEqual(lines[1].strip(), "chr1\t499\t600\tGENE3")

    def test_exon_feature(self):
        """Test extracting exons with a custom name field."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="exon",
            name_key="gene_name",
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0].strip(), "chr1\t99\t150\tGENE1")
        self.assertEqual(lines[1].strip(), "chr1\t149\t200\tGENE1")

    def test_regex_include_patterns(self):
        """Test filtering with regex include patterns."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            filter_key="gene_type",
            bed_columns="detailed",
            include_patterns=[r".*coding$"],  # Matches "protein_coding"
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].strip(), "chr1\t99\t200\t+\tGENE1\tprotein_coding")

    def test_regex_exclude_patterns(self):
        """Test filtering with regex exclude patterns."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            filter_key="gene_type",
            bed_columns="detailed",
            exclude_patterns=[r"pseudo.*"],  # Matches "pseudogene"
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0].strip(), "chr1\t99\t200\t+\tGENE1\tprotein_coding")
        self.assertEqual(lines[1].strip(), "chr1\t299\t400\t+\tGENE2\tlincRNA")

    def test_regex_include_and_exclude_patterns(self):
        """Test filtering with both regex include and exclude patterns."""
        gtf_to_bed(
            input_gtf=self.sample_gtf,
            output_file=self.output_bed,
            feature_type="gene",
            filter_key="gene_type",
            bed_columns="detailed",
            include_patterns=[r".*RNA$"],  # Matches "lincRNA"
            exclude_patterns=[r"pseudo.*"],  # Excludes "pseudogene"
        )
        with open(self.output_bed, "r") as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].strip(), "chr1\t299\t400\t+\tGENE2\tlincRNA")


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.TestLoader().discover("."))
