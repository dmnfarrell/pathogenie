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
import pylab as plt
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from . import tools

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
dbdir = os.path.join(module_path, 'db')

def fetch_sequence_db(name='card'):
    """get sequences"""

    path = dbdir
    links = {'card':'https://github.com/tseemann/abricate/raw/master/db/card/sequences',
            'resfinder':'https://raw.githubusercontent.com/tseemann/abricate/master/db/resfinder/sequences',
            'vfdb':'https://raw.githubusercontent.com/tseemann/abricate/master/db/vfdb/sequences'
            }
    if name in links:
        url = links[name]
    else:
        print('no such name')
        return

    filename = os.path.join(path,"%s.fa" %name)
    if not os.path.exists(filename):
        urllib.request.urlretrieve(url, filename)
    return

def make_blast_database(filenames):
    """Make blast dbs of multiple input files"""

    rec=[]
    for n in filenames:
        seqs = list(SeqIO.parse(n,'fasta'))
        for s in seqs:
            s.id = n + '~' + s.id
        rec.extend(seqs)
    #ref = list(SeqIO.parse('genomes/ecoli_k12.fa','fasta'))
    #ref[0].id = 'ecoli_k12~1'
    #rec.extend(ref)
    SeqIO.write(rec, 'targets.fasta', 'fasta')
    cmd = 'makeblastdb -dbtype nucl -in targets.fasta'
    subprocess.check_output(cmd, shell=True)
    return

def find_genes(target, ref='card', ident=90, coverage=75):
    """Find ref genes by blasting the target sequences"""

    path = os.path.join(dbdir,'%s.fa' %ref)
    dbseqs = list(SeqIO.parse(path,'fasta'))
    print ('blasting %s sequences' %len(dbseqs))
    bl = tools.blast_sequences(target, dbseqs, maxseqs=100, evalue=.1,
                               cmd='blastn', show_cmd=True)
    #print (bl[:5])
    bl['qlength'] = bl.sequence.str.len()
    bl['coverage'] = bl.length/bl.qlength*100
    bl = bl[bl.coverage>coverage]
    bl = bl[bl.pident>ident]
    bl['filename'] = bl.sseqid.apply(lambda x: x.split('~')[0],1)
    bl['id'] = bl.filename.apply(lambda x: os.path.basename(x),1)
    bl['contig'] = bl.sseqid.apply(lambda x: x.split('~')[1],1)
    bl['gene'] = bl['qseqid'].apply(lambda x: x.split('~~~')[1],1)
    
    bl = bl.sort_values('coverage', ascending=False).drop_duplicates(['sstart','send'])
    #print (bl)
    cols = ['qseqid','pident','sstart','send','coverage','contig','gene','id','filename']
    bl = bl[cols]
    return bl

def get_gene_hits(res, gene, filename, db='card'):
    """Get blast hit results"""

    path = os.path.join(dbdir,'%s.fa' %db)
    #dbseqs = SeqIO.to_dict(SeqIO.parse(path,'fasta'))
    dbseqs = tools.fasta_to_dataframe(path)
    dbseqs['gene'] = dbseqs.description.apply(lambda x: x.split('~~~')[1],1)
    #print (dbseqs)
    x = res[res.gene==gene]

    found=[]
    contigs = []
    for i,r in x.iterrows():
        name = r.id
        print (name)
        #if name not in isolates: continue
        seqs = SeqIO.to_dict(SeqIO.parse(r.filename,'fasta'))
        node = r.contig
        if r.sstart<r.send:
            s = seqs[node].seq[r.sstart:r.send]
        else:
            s = seqs[node].seq[r.send:r.sstart].reverse_complement()

        s = SeqRecord(id=name,seq=s)
        found.append(s)
        print (name, r.gene, r['coverage'], r['pident'], len(s), node)
        #add card seq
        contigs.append(seqs[node])

    row = dbseqs[dbseqs.gene==gene].iloc[0]
    print (row)
    found.append(SeqRecord(id=row['name'],seq=Seq(row.sequence)))
    seqfile = 'temp.fa'
    SeqIO.write(found, seqfile,'fasta')
    SeqIO.write(contigs,'contigs.fa','fasta')
    #maaft_alignment(seqfile)
    aln = tools.clustal_alignment(seqfile)
    tools.show_alignment(aln)
    return

def pivot_blast_results(bl):
    x = bl.drop_duplicates(['sstart'])
    m = pd.pivot_table(x, index='id', columns='gene', values='pident')#, aggfunc=np.size)
    #m = m[m.columns[m.loc['ecoli_k12'].isnull()]]
    #m = m.drop('ecoli_k12')
    return m

def plot_heatmap(m, ax=None):

    np.array_split(m, 3)
    if ax == None:
        f,ax=plt.subplots(figsize=(15,8))
    im = ax.imshow(m)
    ax.set_xticks(np.arange(len(m.T)))
    ax.set_yticks(np.arange(len(m)))
    ax.set_xticklabels(m.columns)
    ax.set_yticklabels(m.index)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
    return

def run(filenames=[], db='card', **kwargs):
    """Run pipeline"""

    fetch_sequence_db(db)
    make_blast_database(filenames)
    bl = find_genes('out.fasta', db)
    #find_gene_hits(bl, 'dfrA1_9', '../test_files/RF15B.fa', db)
    bl.to_csv('%s_results.csv' %db)
    m = pivot_blast_results(bl)
    print (m)
    plot_heatmap(m)
    m.to_csv('%s_matrix.csv' %db)
    return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='AMRfinder tool')

    parser.add_argument("-f", "--fasta", dest="filenames", nargs='*',
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-p", "--path", dest="path",
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-d", "--db", dest="db", default='card',
                        help="input fasta file")
    parser.add_argument("-i", "--ident", dest="identity", default='card',
                        help="identity threshold")

    args = vars(parser.parse_args())
    if args['path'] != None:
        args['filenames'] = glob.glob(os.path.join(args['path'],'*.fa*'))
    if len(args['filenames']) == 0:
        return
    run(**args)

if __name__ == '__main__':
    main()
