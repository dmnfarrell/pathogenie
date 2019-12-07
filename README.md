[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<img align="right" src=img/logo.png width=180px>

# pyamrfinder

This is a Python package for detecting antibiotic resistance genes in nucleotide sequences. The inputs are fasta files.

## Usage

pyamrfinder -p <path> -f <filenames> -d <database>

This program utilises the sequence databases compiled by abricate.
    
possible database names:
    
* card
* resfinder
* arg-annot
    
## Links
    
* [abricate by Torsten Seemann](https://github.com/tseemann/abricate)

## Other requirements

* ncbi-blast+ tools 
* clustal