from setuptools import setup
import sys,os

with open('pygenefinder/description.txt') as f:
    long_description = f.read()

setup(
    name = 'pygenefinder',
    version = '0.1.0',
    description = 'find amr genes from fasta sequences',
    long_description = long_description,
    url='https://github.com/dmnfarrell/pygenefinder',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['pygenefinder'],
    package_data={'pygenefinder': ['data/*.*',
                  'description.txt']
                 },
    install_requires=['numpy>=1.10',
                      'pandas>=0.24',
                      'pandastable>=0.12.0',
                      'matplotlib>=3.0',
                      'biopython>=1.5',
                      'future'],
    entry_points = {
        'console_scripts': [
            'pygenefinder=pygenefinder.app:main',
            'pygenefindergui=pygenefinder.gui:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.7',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords = ['bioinformatics','biology','genomics']
)
