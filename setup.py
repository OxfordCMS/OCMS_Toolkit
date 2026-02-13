import sysconfig
import sys
import os
import subprocess
import re
from setuptools import setup, find_packages

setup(
    # package information
    name='ocms_toolkit',
    # NOTE: Keep version format as version="x.y.z" (double quotes, no spaces) for workflow compatibility
    version="1.1.6",
    description='OCMS_Toolkit : Oxford Centre for Microbiome Studies bioinformatic scripts',
    author='Sandi Yen, Nicholas Ilott, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="OCMS, bioinformatics, python",
    url="https://github.com/OxfordCMS/OCMS_Toolkit",
    #packages=find_packages("./") + find_packages("./ocmstoolkit/"),
    packages=find_packages(),
    entry_points={
        'console_scripts': ['ocms_toolkit = ocmstoolkit.ocms_toolkit:main']
    },
    include_package_data=True,
    python_requires='>=3.8.2'
)
