import unittest
import os
from ocmstoolkit.modules.Utility import MetaBam

# global datadir
datadir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
datadir = os.path.join(os.path.dirname(datadir), "fastqs")

# global filename
sam = os.path.join(datadir, "test.sam")
bam = os.path.join(datadir, "test.bam")
cram = os.path.join(datadir, "test.cram")

class TestMetaBam(unittest.TestCase):

    def test_metabam_sam(self):
        mb = MetaBam(sam)
        self.assertEqual(mb.samfile, sam)
        self.assertEqual(mb.bamfile, bam)
        self.assertEqual(mb.cramfile, cram)
        
    def test_metabam_bam(self):
        # Should not raise an exception
        mb = MetaBam(bam)
        self.assertEqual(mb.samfile, sam)
        self.assertEqual(mb.bamfile, bam)
        self.assertEqual(mb.cramfile, cram)

    def test_metabam_cram(self):
        # Should not raise an exception
        mb = MetaBam(cram)
        self.assertEqual(mb.samfile, sam)
        self.assertEqual(mb.bamfile, bam)
        self.assertEqual(mb.cramfile, cram)

if __name__ == '__main__':
    unittest.main()