'''
ocms_toolkit.py - Oxford Centre for Microbiome Studies Toolkit
==============================================================
'''

import os
import sys
import re
import glob
import imp
import cgatcore.iotools as iotools

def main(argv=None):

    argv = sys.argv

    path = os.path.join(os.path.abspath(os.path.dirname(__file__)))
    paths = [path, os.path.abspath(__file__)[:-3]]
    
    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        print((globals()["__doc__"]))

    command = argv[1]
    
    if command in ["new_project", "rename_and_link"]:
        pipeline = command
    else:
        pipeline = "pipeline_{}".format(command)
    
    # remove 'ocms_toolkit' from sys.argv
    del sys.argv[0]

    (file, pathname, description) = imp.find_module(pipeline, paths)
    module = imp.load_module(pipeline, file, pathname, description)
    module.main(sys.argv)
    

if __name__ == "__main__":
    sys.exit(main())
