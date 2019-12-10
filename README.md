[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<img align="right" src=img/logo.png width=180px>

# pyamrfinder

This is a Python package for detecting antibiotic resistance (AMR) genes in nucleotide sequences. It uses Blast to find hits to known AMR gene sequences. The inputs are fasta files.

## Usage

From the command line:

```pyamrfinder -p <path-to-fasta-files> -d <database>```

or 

```pyamrfinder -f <filename> -d <database>```

This program utilises the sequence databases compiled by abricate. Possible database names:
    
* card
* resfinder
* arg-annot
* resfinder
* ncbi
* ecoh

You can also use a graphical application. There is an installer provided for windows users.

<img src=img/screenshot1.png width=480px>

## Installation

Any operating systems:

```pip install pyamrfinder```

You also require ncbi-blast+ tools. This can be installed on Debian/Ubuntu based systems using:

```sudo apt install ncbi-blast+```

On windows you can install the binary from [this link](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-win64.exe).

### Windows

Most windows users will want to use the installer for the graphical application.

## Links
    
* [abricate by Torsten Seemann](https://github.com/tseemann/abricate)
