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

You can also use a graphical application. It can be launched from the terminal using:

```pyamrfindergui```

<img src=img/screenshot1.png width=480px>

There is a self-contained graphical application for windows users (see below).

## Installation

All operating systems with Python (>=3.6 recommended) installed:

```pip install pyamrfinder```

### Windows GUI

Most windows users will probably want to use the bundled graphical application. [Download here](https://github.com/dmnfarrell/pyamrfinder/releases/download/v0.1.0/pyamrfinder-0.1.0-win64.zip). Just unzip to a folder and launch the program **pyamrfinder.exe**.

This has all the dependencies bundled with the program.

## Dependencies

You also require ncbi-blast+ tools and clustalw. These can be installed on Debian/Ubuntu based systems using:

```sudo apt install ncbi-blast+ clustal```

On Windows, (if installing via pip) you can install blast tools from here:

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-win64.exe

## Links

* [abricate by Torsten Seemann](https://github.com/tseemann/abricate)
