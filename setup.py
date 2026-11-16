from setuptools import setup, find_packages

setup(
    # package information
    name='ocms_toolkit',
    # NOTE: Keep version format as version="x.y.z" (double quotes, no spaces) for workflow compatibility
    version="1.1.7",
    description='OCMS_Toolkit : Oxford Centre for Microbiome Studies bioinformatic scripts',
    author='Sandi Yen, Nicholas Ilott, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="OCMS, bioinformatics, python",
    url="https://github.com/OxfordCMS/OCMS_Toolkit",
    packages=find_packages(),
    # Snakemake workflows and their configuration templates
    package_data={
        "ocmstoolkit": ["workflows/*.smk", "config/*.yaml"],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': ['ocms_toolkit = ocmstoolkit.cli:main']
    },
    install_requires=[
        "pysam",
        "pyyaml",
        "snakemake>=9",
        "snakemake-executor-plugin-slurm",
    ],
    extras_require={
        "test": ["pytest"],
    },
    # Snakemake 9 requires Python 3.11 or later
    python_requires='>=3.11'
)
