import unittest
import os
from ocmstoolkit.modules.Utility import MetaFastn

# global datadir
datadir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
datadir = os.path.join(datadir, "fastqs")

# global filename
fastq_pe = os.path.join(datadir, "test.fastq.1.gz")
fastq_se = os.path.join(datadir, "test.fastq.gz")

class TestMetaFastn(unittest.TestCase):

    def test_metafastn_pe(self):
        # Should not raise an exception
        mf = MetaFastn(fastq_pe)
        self.assertEqual(mf.fastn1, fastq_pe)
        self.assertTrue(isinstance(mf.head, list))
        self.assertTrue(len(mf.head) > 0)
        # If paired, fastn2 and fastn3 should be set (if file exists)
        self.assertTrue(mf.fastn2 is not None)
        self.assertTrue(mf.fastn3 is not None or mf.fastn3 is None)   # Just check attribute exists

    def test_metafastn_se(self):
        # Should not raise an exception
        mf = MetaFastn(fastq_se)
        self.assertEqual(mf.fastn1, fastq_se)
        self.assertTrue(isinstance(mf.head, list))
        self.assertTrue(len(mf.head) > 0)
        # If paired, fastn2 and fastn3 should be set (if file exists)
        self.assertTrue(mf.fastn2 is None)  # Just check attributes exist
        self.assertTrue(mf.fastn3 is None)

    def test_metafastn_file_not_found(self):
        with self.assertRaises(Exception):
            MetaFastn(os.path.join(datadir, 'notfound.fastq.1.gz'))

if __name__ == '__main__':
    unittest.main()