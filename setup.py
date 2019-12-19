import os
import io
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def get_version(relpath):
    '''Read version info from a file without importing it'''
    for line in io.open(os.path.join(os.path.dirname(__file__), relpath), encoding='cp437'):
        if '__version__' in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]

setup(
    name = "sag_mg_recruit",
    version = get_version('sag_mg_recruit.py'),
    author = "Julia Brown",
    author_email = "julia@bigelow.org",
    description = ("recruit metagenomes to SAGs"),
    py_modules=['sag_mg_recruit'],
    license = "MIT",
    install_requires=[
        'click',
        'pandas',
        'pysam',
        'sarge',
        'matplotlib',
        'biopython'],
    entry_points='''
        [console_scripts]
        sag-mg-recruit=sag_mg_recruit:main
        '''
)
