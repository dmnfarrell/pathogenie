#!/usr/bin/env python

"""
    AMR finder.
    Created Nov 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

from __future__ import absolute_import, print_function
import sys,os,subprocess,glob
import urllib
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from . import tools

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')

def fetch_sequence_db(name='card'):
    """get sequences"""

    path = datadir
    if name == 'card':
        url = 'https://github.com/tseemann/abricate/raw/master/db/card/sequences'
    filename = os.path.join(path,"%s.fa" %name)
    if not os.path.exists(filename):
        urllib.request.urlretrieve(url, filename)
    return

def blast_card(ident=90,coverage=.75):
    """blast card seqs"""

    cardseqs = list(SeqIO.parse('/local/abricate/db/card/sequences','fasta'))
    #bl = tools.local_blast(qfile, 'scaffolds.fasta', ident=ident, params='-e 10 -a 2 -S 1')
    bl = tools.blast_sequences('scaffolds.fasta',cardseqs,maxseqs=100,evalue=.1,cmd='blastn',show_cmd=True)
    print (bl[:5])
    bl['qlength'] = bl.sequence.str.len()
    bl['coverage'] = bl.length/bl.qlength
    bl = bl[bl.coverage>coverage]
    bl = bl[bl.pident>ident]
    bl['id'] = bl.sseqid.apply(lambda x: x.split('~')[0],1)
    bl['contig'] = bl.sseqid.apply(lambda x: x.split('~')[1],1)
    bl['gene'] = bl['qseqid'].apply(lambda x: x.split('~~~')[1],1)
    #print (bl)
    cols = ['qseqid','sseqid','pident','sstart','send','coverage','contig','gene','id']
    bl = bl[cols]
    return bl

def make_blast_database(filenames):
    """make blast dbs"""

    rec=[]
    for n in pig:
        seqs = list(SeqIO.parse('scaffolds/%s.fa' %n,'fasta'))
        for s in seqs:
            s.id = n + '~' + s.id
        rec.extend(seqs)
    ref = list(SeqIO.parse('genomes/ecoli_k12.fa','fasta'))
    ref[0].id = 'ecoli_k12~1'
    rec.extend(ref)
    SeqIO.write(rec, 'scaffolds.fasta', 'fasta')
    cmd = 'makeblastdb -dbtype nucl -in scaffolds.fasta'

def find_hits(res, gene, db='card'):
    """get blast hit results"""

    x = res[res.gene==gene]
    found=[]
    nodes=[]
    for i,r in x.iterrows():
        name=r.id
        if name not in isolates: continue
        seqs = SeqIO.to_dict(SeqIO.parse('scaffolds/%s.fa' %r.id,'fasta'))
        node = r.contig
        if r.sstart<r.send:
            s = seqs[node].seq[r.sstart:r.send]
        else:
            s = seqs[node].seq[r.send:r.sstart].reverse_complement()

        s = SeqRecord(id=name,seq=s)
        found.append(s)
        print (name, r.gene, r['coverage'], r['pident'], len(s), node)
        #add card seq
        nodes.append(seqs[node])

    found.append(cardseqs[gene])
    seqfile = 'card_%s.fa' %gene
    SeqIO.write(found,seqfile,'fasta')
    SeqIO.write(nodes,'ctx-contigs.fa','fasta')
    #maaft_alignment(seqfile)
    tools.clustal_alignment(seqfile)

def run(filename):
    fetch_sequence_db()

    return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='AMRfinder tool')

    parser.add_argument("-f", "--fasta", dest="filename",
                        help="input fasta file", metavar="FILE")
    args = vars(parser.parse_args())
    run(**args)

if __name__ == '__main__':
    main()
