from unittest import TestCase
import ocmstoolkit.modules.Utility as Utility
import os

# global datadir
datadir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
datadir = os.path.join(datadir, "fastqs")

# global filename
fastq = os.path.join(datadir, "test.fastq.1.gz")

# Note that in the test we are using the infile and
# outfile as the same file
class TestBaseTool(TestCase):

    def test_base_tool_infile(self):
        fastq_btool = Utility.BaseTool(fastq, fastq)
        self.assertEqual(fastq_btool.infile, fastq)

    def test_base_tool_outfile(self):
        fastq_btool = Utility.BaseTool(fastq, fastq)
        self.assertEqual(fastq_btool.outfile, fastq)

    def test_base_tool_indir(self):
        indir = datadir
        fastq_btool = Utility.BaseTool(fastq, fastq)
        self.assertEqual(fastq_btool.indir, indir)

    def test_base_tool_outdir(self):
        outdir = datadir
        fastq_btool = Utility.BaseTool(fastq, fastq)
        self.assertEqual(fastq_btool.outdir, outdir)

    def test_base_tool_params(self):
        fastq_btool = Utility.BaseTool(fastq, fastq, attribute="fastq")
        self.assertEqual(fastq_btool.PARAMS, {"attribute": "fastq"})
        
