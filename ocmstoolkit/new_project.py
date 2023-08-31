'''
new_project.py
==============

:Author: Sandi Yen
:Tags: Python

Overview
========

Use this script to initiate a new project. New directories and symlinked directories will be created

Usage
-----

The script takes two parameters. First parameter `project_name` is the name of the new project. Second paramter `level` is the level for which project should be set up: group, user, or both. Setting `level=group` sets up directories necessary in the project folder. Setting `level=user` sets up directories necessary in the user directory. Setting `level=both` sets up all directories and symlinks.

Example::
    ocms new_project --project_name=<PROJECT-NAME> --level=both
    ocms new_project -p <PROJECT-NAME>
    ocms new_project -p <PROJECT-NAME> -l user
    
Type::
  
    ocms new_project --help

for comand line help.

'''

import sys
import os
import getpass
import pwd
import grp
import cgatcore.experiment as E

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """
    
    if argv is None:
        argv = sys.argv
    
    # set up command line parser
    parser = E.ArgumentParser(description=__doc__)
    
    parser.add_argument("-p", "--project_name", dest="project_name", type=str,
                        required=True, help="name of project")
    parser.add_argument("-l", "--level", dest="level", type=str,
                        choices=["group","user","both"], default="both",
                        help="level at which to set up project")
    
    # add common options (-h/--help) and parse command line
    (args) = E.start(parser, argv=argv)
    
    # get local username and group
    username = getpass.getuser()
    groups = [g.gr_name for g in grp.getgrall() if username in g.gr_mem]
    gid = pwd.getpwnam(username).pw_gid
    groups.append(grp.getgrgid(gid).gr_name)
    # just use first group
    group = groups[0]
    
    # set up output directories
    project_dir = "/gpfs3/well/" + group + "/projects/" + args.project_name
    archive_dir = "/gpfs3/well/" + group + "/projects/archive/" + args.project_name
    data_dir = "/gpfs3/well/" + group + "/projects/" + args.project_name + "/data"
    work_dir = "/gpfs3/well/" + group + "/users/" + username + "/work/" + args.project_name
    devel_dir = "/gpfs3/well/" + group + "/users/" + username + "/devel/" + args.project_name
    code_dir = "/gpfs3/well/" + group + "/users/" + username + "/devel/" + args.project_name + "/code"
    analysis_dir = "/gpfs3/well/" + group + "/users/" + username + "/work/" + args.project_name + "/analysis"
    
    
    # make directories, failing if output directories already exist
    # project directories
    if args.level in ['group','both']:
    
        os.makedirs(project_dir, exist_ok=False)
        os.makedirs(data_dir, exist_ok=False)
        os.makedirs(archive_dir, exist_ok=False)
        
    # work and devel directories
    if args.level in ['user','both']:
        os.makedirs(work_dir, exist_ok=False)
        os.makedirs(analysis_dir, exist_ok=False)
        os.makedirs(devel_dir, exist_ok=False)
        os.makedirs(code_dir, exist_ok=False)
        
        # symlink code
        os.symlink(code_dir, work_dir + "/code")
            
        # write footer and output benchmark information
        E.stop()
        
if __name__ == "__main__":
    sys.exit(main(sys.argv))