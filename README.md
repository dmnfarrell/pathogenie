[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# pathogenie

<img align="right" src=img/logo.png width=180px>

This is a desktop and command line program for annotating draft bacterial and viral genomes. It may also be used for quickly detecting arbitrary sequences such as antibiotic resistance genes (AMR) proteins in nucleotide sequences. It uses Blast to find hits to known gene sequences from sequence databases. The inputs are fasta files. Annotation is performed in a similar manner to Prokka and first  requires an assembled genome if you have sequenced reads. The program is written in Python.

## Usage

From the command line:

```pathogenie -p <path-to-fasta-files> -d <database>```

or

```pathogenie -f <filename> -d <database>```

This program utilises the sequence databases compiled by abricate. Possible database names:

* card
* resfinder
* arg-annot
* resfinder
* ncbi
* ecoh

You can also use a graphical application. It can be launched from the terminal using:

```pathogenie-gui```

<img src=img/screenshot1.png width=480px>

There is a self-contained graphical application for windows users (see below).

## Installation

All operating systems with Python (>=3.6 required) installed:

```pip install -e git+https://github.com/dmnfarrell/pathogenie.git#egg=pathogenie```

## Dependencies

You require ncbi-blast+ tools and clustalw for basic gene finding. The following programs are used for genome annotation:

* prodigal
* hmmer3
* aragorn

Windows: these will be downloaded for you when you launch the program.

Linux: these can all be installed on Debian/Ubuntu based systems using:

```apt install ncbi-blast+ clustal prodigal aragorn hmmer```

OSX: **Not tested** but may work if you can install the dependencies. You can probably install them with bioconda.

### Windows GUI

Most windows users will probably want to use the bundled graphical application. [Download current version here](https://github.com/dmnfarrell/pathogenie/releases/download/0.1.0/pathogenie-0.1.0-win64.zip). Just unzip to a folder and launch the program **pathogenie.exe**.

This has all the dependencies bundled with the program.

## Links

* [Prokka](https://github.com/tseemann/prokka)
* [abricate](https://github.com/tseemann/abricate)
