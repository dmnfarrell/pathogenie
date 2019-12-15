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
import sys,os,subprocess,glob,re
import urllib
import tempfile
import pandas as pd
import numpy as np
import pylab as plt
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from . import tools

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config/pyamrfinder')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
dbdir = os.path.join(config_path, 'db')
prokkadbdir = os.path.join(config_path, 'prokka')
prokka_db_names = ['sprot','IS','AMR']
links = {'card':'https://github.com/tseemann/abricate/raw/master/db/card/sequences',
        'resfinder':'https://raw.githubusercontent.com/tseemann/abricate/master/db/resfinder/sequences',
        'vfdb':'https://raw.githubusercontent.com/tseemann/abricate/master/db/vfdb/sequences',
        'sprot':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/sprot'}

if not os.path.exists(config_path):
    try:
        os.makedirs(config_path, exist_ok=True)
    except:
        os.makedirs(config_path)

db_names = ['card','resfinder','argannot','ncbi','plasmidfinder','ecoh','vfdb']

def check_databases():
    """download databases"""

    print ('checking databases')
    os.makedirs(dbdir, exist_ok=True)
    for name in db_names:
        #print (name)
        fetch_sequence_from_url(name)
    return

def fetch_sequence_from_url(name='card', path=None):
    """get sequences"""

    if path == None:
        path = dbdir
    if not os.path.exists(path):
        os.makedirs(path)

    if name in links:
        url = links[name]
    else:
        print('no such name')
        return

    filename = os.path.join(path,"%s.fa" %name)
    if not os.path.exists(filename):
        urllib.request.urlretrieve(url, filename)
    return

def make_target_database(filenames):
    """Make blast db from multiple input files"""

    rec=[]
    for n in filenames:
        seqs = list(SeqIO.parse(n,'fasta'))
        for s in seqs:
            s.id = n + '~' + s.id
        rec.extend(seqs)

    targfile = os.path.join(tempdir, 'targets.fasta')
    SeqIO.write(rec, targfile, 'fasta')
    tools.make_blast_database(targfile)
    return

def find_genes(target, ref='card', ident=90, coverage=75, duplicates=False, **kwds):
    """Find ref genes by blasting the target sequences"""

    path = os.path.join(dbdir,'%s.fa' %ref)
    #the AMR db is the query for the blast
    queryseqs = list(SeqIO.parse(path,'fasta'))
    print ('blasting %s sequences' %len(queryseqs))
    bl = tools.blast_sequences(target, queryseqs, maxseqs=100, evalue=.1,
                               cmd='blastn', show_cmd=True)

    bl['qlength'] = bl.sequence.str.len()
    bl['coverage'] = bl.length/bl.qlength*100
    bl = bl[bl.coverage>coverage]
    bl = bl[bl.pident>ident]
    bl['filename'] = bl.sseqid.apply(lambda x: x.split('~')[0],1)
    bl['id'] = bl.filename.apply(lambda x: os.path.basename(x),1)
    bl['contig'] = bl.sseqid.apply(lambda x: x.split('~')[1],1)
    bl['gene'] = bl['qseqid'].apply(lambda x: x.split('~~~')[1],1)

    #remove exact and close duplicates
    bl = bl.sort_values(['coverage','pident'], ascending=False).drop_duplicates(['contig','sstart','send'])
    if duplicates == False:
        dist = 20
        bl=bl.sort_values(by=["contig","sstart"])
        unique = bl.sstart.diff().fillna(dist)
        bl = bl[unique>=dist]
    cols = ['gene','id','qseqid','pident','coverage','sstart','send','contig','filename']
    bl = bl[cols]
    return bl

def get_gene_hits(res, gene, filename, db='card'):
    """Get blast hit results for a gene"""

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
        #print (name)
        #if name not in isolates: continue
        seqs = SeqIO.to_dict(SeqIO.parse(r.filename,'fasta'))
        node = r.contig
        if r.sstart<r.send:
            s = seqs[node].seq[r.sstart:r.send]
        else:
            s = seqs[node].seq[r.send:r.sstart].reverse_complement()

        s = SeqRecord(id=name,seq=s)
        found.append(s)
        #print (name, r.gene, r['coverage'], r['pident'], len(s), node)
        #add card seq
        contigs.append(seqs[node])

    row = dbseqs[dbseqs.gene==gene].iloc[0]
    found.append(SeqRecord(id=row['name'],seq=Seq(row.sequence)))
    SeqIO.write(contigs,'contigs.fa','fasta')
    return found

def get_alignment(seqs):

    seqfile = 'temp.fa'
    SeqIO.write(seqs, seqfile,'fasta')
    aln = tools.clustal_alignment(seqfile)
    tools.show_alignment(aln)
    return aln

def pivot_blast_results(bl):

    x = bl.drop_duplicates(['sstart'])
    m = pd.pivot_table(x, index='id', columns='gene', values='pident')#, aggfunc=np.size)
    #m = m[m.columns[m.loc['ecoli_k12'].isnull()]]
    #m = m.drop('ecoli_k12')
    return m

def plot_heatmap(m, fig=None, title=''):

    from matplotlib.gridspec import GridSpec
    l=1+int(len(m)/30)
    if fig == None:
        fig = plt.figure()
    gs = fig.add_gridspec(1, l)
    chunks = np.array_split(m,l)
    i=0
    for df in chunks:
        ax = fig.add_subplot(gs[0,i])
        im = ax.imshow(df, cmap='Blues')
        ax.set_xticks(np.arange(len(df.T)))
        ax.set_yticks(np.arange(len(df)))
        ax.set_xticklabels(df.columns)
        ax.set_yticklabels(df.index)
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
                 rotation_mode="anchor")
        i+=1
    fig.suptitle(title)
    fig.subplots_adjust(hspace=1.2, bottom=.2)
    return

def prodigal(infile):
    """Run prodigal"""

    cmd = 'prodigal'
    if getattr(sys, 'frozen', False):
        cmd = tools.resource_path('bin/prodigal.exe')
    name = os.path.splitext(infile)[0]
    cmd = '{c} -i {n}.fa -a {n}.faa -f gff -o {n}.gff -p single'.format(c=cmd,n=name)
    subprocess.check_output(cmd, shell=True)
    resfile = name+'.faa'
    return resfile

def get_prodigal_coords(x):
    s = re.split('\#|\s',x.replace(' ',''))
    coords = [int(i) for i in s[1:4]]
    return  pd.Series(coords)

def prokka_header_info(x):
    s = re.split('~~~',x)
    return pd.Series(s)

def annotate_contigs(infile, outfile=None, **kwargs):
    """
    Annotate nucelotide sequences (usually a draft assembly with contigs)
    using prodigal and blast to prokka seqs. Writes a genbank file to the
    same folder.
    Args:
        infile: input fasta file
        outfile: output genbank
    returns:
        a list of SeqRecords with the features
    """

    #run prodigal
    resfile = prodigal(infile)
    print (resfile)
    #get target seqs
    seqs = list(SeqIO.parse(resfile,'fasta'))
    #make blast db of prokka proteins
    dbname = os.path.join(prokkadbdir,'sprot.fa')
    tools.make_blast_database(dbname, dbtype='prot')
    print ('blasting ORFS to uniprot sequences')
    bl = tools.blast_sequences(dbname, seqs, maxseqs=100, evalue=.01,
                                cmd='blastp', show_cmd=True, **kwargs)

    bl[['protein_id','gene','product','cog']] = bl.stitle.apply(prokka_header_info,1)

    cols = ['qseqid','sseqid','pident','sstart','send','protein_id','gene','product']
    x = bl.sort_values(['qseqid','pident'], ascending=False).drop_duplicates(['qseqid'])[cols]

    #read input file seqs
    contigs = SeqIO.to_dict(SeqIO.parse(infile,'fasta'))
    #read in prodigal fasta to dataframe
    df = tools.fasta_to_dataframe(resfile)

    df[['start','end','strand']] = df.description.apply(get_prodigal_coords,1)
    #merge blast result with prodigal fasta file info
    df = df.merge(x,left_on='name',right_on='qseqid',how='left')
    df['contig'] = df['name'].apply(lambda x: x[:6])

    def get_contig(x):
        return ('_').join(x.split('_')[:-1])
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Alphabet import generic_dna
    l=1
    if outfile is None:
        outfile = infile+'.gbk'
    handle = open(outfile,'w+')
    recs = []
    #group by contig and get features for each protein found
    for c,df in df.groupby('contig'):
        contig = get_contig(df.iloc[0]['name'])
        nucseq = contigs[contig].seq
        rec = SeqRecord(nucseq,id=c)
        rec.seq.alphabet = generic_dna
        for i,row in df.iterrows():
            tag = 'PREF_{l:04d}'.format(l=l)
            quals = {'gene':row.gene,'product':row['product'],'locus_tag':tag,'translation':row.sequence}
            feat = SeqFeature(FeatureLocation(row.start,row.end, row.strand), strand=row.strand,
                              type="CDS", qualifiers=quals)
            rec.features.append(feat)
            l+=1
        #print(rec.format("gb"))
        SeqIO.write(rec, handle, "genbank")
        recs.append(rec)
    handle.close()
    return recs

def run(filenames=[], db='card', outdir='amr_results', **kwargs):
    """Run pipeline"""

    check_databases()
    make_target_database(filenames)
    targfile = os.path.join(tempdir, 'targets.fasta')
    print (targfile)
    bl = find_genes(targfile, db, **kwargs)
    m = pivot_blast_results(bl)
    #plot_heatmap(m)
    if outdir != None:
        os.makedirs(outdir, exist_ok=True)
        bl.to_csv(os.path.join(outdir,'%s_results.csv' %db))
        m.to_csv(os.path.join(outdir,'%s_matrix.csv' %db))
    print ('results saved to %s' %outdir)
    return

def run_test():
    files = glob.glob(os.path.join(datadir, '*.fa'))
    run(files)
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
    parser.add_argument("-o", "--out", dest="outdir", default='amr_results',
                        help="output folder", metavar="FILE")
    parser.add_argument("-d", "--db", dest="db", default='card',
                        help="input fasta file")
    parser.add_argument("-i", "--ident", dest="identity", default='card',
                        help="identity threshold")
    parser.add_argument("-t", "--test", dest="test", action='store_true',
                        help="test run")

    args = vars(parser.parse_args())
    if args['test'] == True:
        run_test()
        return
    if args['path'] != None:
        args['filenames'] = glob.glob(os.path.join(args['path'],'*.fa*'))
    if args['filenames'] == None:
        return
    run(**args)

if __name__ == '__main__':
    main()
