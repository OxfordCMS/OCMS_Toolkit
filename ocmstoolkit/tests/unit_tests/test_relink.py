import unittest
import os
import tempfile
from ocmstoolkit.modules.Utility import relink

class TestRelink(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.src_file = os.path.join(self.temp_dir.name, 'source.txt')
        self.dst_link = os.path.join(self.temp_dir.name, 'link.txt')
        with open(self.src_file, 'w') as f:
            f.write("test")

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_symlink_creation(self):
        relink(self.src_file, self.dst_link)
        self.assertTrue(os.path.islink(self.dst_link))
        self.assertEqual(os.readlink(self.dst_link), os.path.abspath(self.src_file))

    def test_replaces_existing_link(self):
        # First create the link
        relink(self.src_file, self.dst_link)
        # Now, create a dummy file at the link path
        with open(self.dst_link, 'w') as f:
            f.write("should be replaced")
        # Relink should replace file with symlink
        relink(self.src_file, self.dst_link)
        self.assertTrue(os.path.islink(self.dst_link))
        self.assertEqual(os.readlink(self.dst_link), os.path.abspath(self.src_file))

if __name__ == '__main__':
    unittest.main()
