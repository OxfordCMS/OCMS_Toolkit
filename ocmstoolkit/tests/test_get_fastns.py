from tempfile import NamedTemporaryFile
from unittest import TestCase
import ocmstoolkit.modules.Utility as Utility
import os

# global datadir
datadir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
datadir = os.path.join(datadir, "fastqs")

class TestGetFastns(TestCase):

    def test_fastn1s(self):
        fastns = Utility.get_fastns(datadir, 1)
        self.assertEqual(fastns[0], os.path.join(datadir, "test.fastq.1.gz"))

    def test_fastn2s(self):
        fastns = Utility.get_fastns(datadir, 1,2)
        self.assertEqual(fastns, tuple([os.path.join(datadir, x)]
                                       for x in ["test.fastq.1.gz", "test.fastq.2.gz"]))

    def test_fastn3s(self):
        fastns = Utility.get_fastns(datadir, 1,2,3)
        self.assertEqual(fastns, tuple([os.path.join(datadir, x)]
                                       for x in ["test.fastq.1.gz",
                                                 "test.fastq.2.gz",
                                                 "test.fastq.3.gz"]))
        
