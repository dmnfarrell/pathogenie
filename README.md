[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# pathogenie

<img align="right" src=img/logo.png width=180px>

This is a desktop and command line program for annotating draft bacterial and viral genomes. It may also be used for quickly detecting arbitrary sequences such as antibiotic resistance genes (AMR) proteins in nucleotide sequences. It uses Blast to find hits to known gene sequences from sequence databases. The inputs are fasta files. Annotation is performed in a similar manner to Prokka and first requires an assembled genome if you have sequenced reads. The program is written in Python. Currently it is available as a graphical desktop application. A command line tool will also be added. You can also used from inside Python.

## Usage

The graphical application can be launched from the terminal using:

```pathogenie-gui```

From the GUI you may load fasta files into a table and then run genome annotation or gene finding with custom databases. This program utilises the sequence databases for gene finding compiled by abricate:

* card
* resfinder
* arg-annot
* resfinder
* ncbi
* ecoh

The GUI layout is shown below:

<img src=img/screenshot1.png width=480px>

## Using in python

Run an annotation on a fasta file like a set of contigs:

```python
import pathogenie
featdf,recs = pathogenie.app.run_annotation(filename, threads=10, kingdom='bacteria')
#save to genbank
pathogenie.tools.recs_to_genbank(recs, gbfile)
```

## Installation

All operating systems with Python (>=3.6 required) installed:

```pip install -e git+https://github.com/dmnfarrell/pathogenie.git#egg=pathogenie```

## Dependencies

You require ncbi-blast+ tools and clustalw for basic gene finding. The following programs are used for genome annotation:

* prodigal
* hmmer3
* aragorn

### Windows

These executables will be downloaded for you when you first launch the program.

Blast also requires the Visual Studio 2015 C++ redistributable runtime package: https://www.microsoft.com/en-us/download/details.aspx?id=48145

### Linux

The external binaries can all be installed on Debian/Ubuntu based systems using:

```sudo apt install ncbi-blast+ clustal prodigal aragorn hmmer```

### OSX

**Not yet tested** but may work if you can install the dependencies. You can probably install them with bioconda.

## Links

* [Prokka](https://github.com/tseemann/prokka)
* [abricate](https://github.com/tseemann/abricate)
