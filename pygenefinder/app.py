#!/usr/bin/env python

"""
    pygenefinder cmd line tool.
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
import urllib, hashlib
import tempfile
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from . import tools

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config/pygenefinder')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
dbdir = os.path.join(config_path, 'db')
customdbdir = os.path.join(config_path, 'custom')
hmmdir = os.path.join(config_path, 'hmms')
prokkadbdir = os.path.join(config_path, 'prokka')
prokka_db_names = ['sprot','IS','AMR']
links = {'card':'https://github.com/tseemann/abricate/raw/master/db/card/sequences',
        'resfinder':'https://raw.githubusercontent.com/tseemann/abricate/master/db/resfinder/sequences',
        'vfdb':'https://raw.githubusercontent.com/tseemann/abricate/master/db/vfdb/sequences',
        'ncbi':'https://raw.githubusercontent.com/tseemann/abricate/master/db/ncbi/sequences',
        'argannot':'https://raw.githubusercontent.com/tseemann/abricate/master/db/argannot/sequences',
        'ecoh':'https://raw.githubusercontent.com/tseemann/abricate/master/db/ecoh/sequences',
        'plasmidfinder':'https://raw.githubusercontent.com/tseemann/abricate/master/db/plasmidfinder/sequences',
        'sprot':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/sprot',
        'amr':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/AMR',
        'IS':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/IS',
        'bacteria.16SrRNA': 'https://raw.githubusercontent.com/dmnfarrell/pygenefinder/master/db/bacteria.16SrRNA.fna',
        'bacteria.23SrRNA': 'https://raw.githubusercontent.com/dmnfarrell/pygenefinder/master/db/bacteria.23SrRNA.fna',
        'HAMAP.hmm':'https://github.com/tseemann/prokka/raw/master/db/hmm/HAMAP.hmm'  }

hmm = 'https://github.com/tseemann/prokka/raw/master/db/hmm/HAMAP.hmm'

if not os.path.exists(config_path):
    try:
        os.makedirs(config_path, exist_ok=True)
    except:
        os.makedirs(config_path)

db_names = ['card','resfinder','argannot','ncbi','plasmidfinder','ecoh','vfdb',
            'bacteria.16SrRNA','bacteria.23SrRNA']

def check_databases():
    """download databases"""

    print ('checking databases')
    os.makedirs(dbdir, exist_ok=True)
    for name in db_names:
        fetch_sequence_from_url(name)
    return

def fetch_sequence_from_url(name='card', path=None, ext='.fa'):
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

    filename = os.path.join(path,"%s%s" %(name,ext))
    print (filename)
    if not os.path.exists(filename):
        try:
            urllib.request.urlretrieve(url, filename)
        except:
            print ('no such URL?')
    return

def make_target_database(filenames):
    """Make blast db from multiple input files"""

    rec=[]
    if filenames is str:
        filenames = [filenames]
    for n in filenames:
        seqs = list(SeqIO.parse(n,'fasta'))
        for s in seqs:
            s.id = n + '~' + s.id
        rec.extend(seqs)

    targfile = os.path.join(tempdir, 'targets.fasta')
    SeqIO.write(rec, targfile, 'fasta')
    tools.make_blast_database(targfile)
    return targfile

def find_genes(target, ref='card', ident=90, coverage=75, duplicates=False, threads=2, **kwds):
    """Find ref genes by blasting the target sequences"""

    path = os.path.join(dbdir,'%s.fa' %ref)
    #the AMR db is the query for the blast
    queryseqs = list(SeqIO.parse(path,'fasta'))
    print ('blasting %s sequences' %len(queryseqs))
    bl = tools.blast_sequences(target, queryseqs, maxseqs=1, evalue=1e-4,
                               cmd='blastn', show_cmd=True, threads=int(threads))

    bl['qlength'] = bl.sequence.str.len()
    bl['coverage'] = bl.length/bl.qlength*100
    bl = bl[bl.coverage>coverage]
    bl = bl[bl.pident>ident]
    bl['filename'] = bl.sseqid.apply(lambda x: x.split('~')[0],1)
    bl['id'] = bl.filename.apply(lambda x: os.path.basename(x),1)
    bl['contig'] = bl.sseqid.apply(lambda x: x.split('~')[1],1)
    try:
        bl['gene'] = bl['qseqid'].apply(lambda x: x.split('~~~')[1],1)
    except:
        bl['gene'] = bl.qseqid

    #remove exact and close duplicates
    print (len(bl))
    bl = bl.sort_values(['bitscore'], ascending=False).drop_duplicates(['contig','sstart','send'])
    print (len(bl))
    if duplicates == False:
        dist = 20
        x=bl.sort_values(by=["contig","sstart"],ascending=False)
        #print (x[:15][x.columns[:5]])
        unique = x.sstart.diff().abs().fillna(dist)
        bl = bl[unique>=dist]
    cols = ['gene','id','qseqid','pident','coverage','sstart','send','contig','description','filename','bitscore']
    bl = bl[cols]
    return bl

def get_gene_hits(res, gene, filename, db='card'):
    """Get blast hit results for a gene"""

    path = os.path.join(dbdir,'%s.fa' %db)
    #dbseqs = SeqIO.to_dict(SeqIO.parse(path,'fasta'))
    dbseqs = tools.fasta_to_dataframe(path)
    try:
        dbseqs['gene'] = dbseqs.description.apply(lambda x: x.split('~~~')[1],1)
    except:
        print (dbseqs)
        dbseqs['gene'] = dbseqs['name']
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

    import pylab as plt
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

def create_locus_tag(filename):

    x = hashlib.md5(filename.encode('utf-8')).hexdigest()
    x = ''.join(i for i in x if not i.isdigit())
    x = x.upper()[:7]
    return x

def prodigal(infile):
    """Run prodigal"""

    cmd = 'prodigal'
    if getattr(sys, 'frozen', False):
        cmd = tools.resource_path('bin/prodigal.exe')
    #name = os.path.splitext(infile)[0]
    name = os.path.join(tempdir,'prodigal')
    cmd = '{c} -i {i} -a {n}.faa -f gff -o {n}.gff -p single'.format(i=infile,c=cmd,n=name)
    subprocess.check_output(cmd, shell=True)
    resfile = os.path.join(tempdir, 'prodigal.faa')
    return resfile

def get_prodigal_coords(x):
    s = re.split('\#|\s',x.replace(' ',''))
    coords = [int(i) for i in s[1:4]]
    return  pd.Series(coords)

def prokka_header_info(x):
    s = re.split('~~~',x)
    return pd.Series(s)

def hmmer(infile, threads=4, hmm_file=None):
    """Run hmmer"""

    def get_contig(x):
        return ('_').join(x.split('_')[:-1])
    if getattr(sys, 'frozen', False):
        hmmpresscmd = tools.resource_path('bin/hmmpress.exe')
        hmmscancmd = tools.resource_path('bin/hmmscan.exe')
    else:
        hmmpresscmd = 'hmmpress'
        hmmscancmd = 'hmmscan'
    df = tools.fasta_to_dataframe(infile)
    fetch_sequence_from_url('HAMAP.hmm', hmmdir, ext='')
    db = os.path.join(hmmdir,'HAMAP.hmm')
    cmd = '%s -f %s' %(hmmpresscmd,db)
    tmp = subprocess.check_output(cmd, shell=True)
    out='hmm.txt'
    cmd = "{c} --noali --notextw --acc -E 1e-4 --cpu {t} --tblout {o} -o /tmp/hmm.out {db} {i}".format(c=hmmscancmd,t=threads,db=db,i=infile,o=out)
    print (cmd)
    tmp = subprocess.check_output(cmd, shell=True)
    h = tools.read_hmmer3('hmm.txt')
    #print (df)
    df = h.merge(df,on='name',how='left')
    #print (df)
    #get coords from prodigal fasta heading if available
    df[['start','end','strand']] = df.description.apply(get_prodigal_coords,1)
    df['contig'] = df['name'].apply(get_contig)
    return df

def aragorn(infile):
    """Run aragorn"""

    cmd = 'aragorn'
    if getattr(sys, 'frozen', False):
        cmd = tools.resource_path('bin/aragorn.exe')
    cmd = '{c} -l -gcbact -t -w {i} -o /tmp/aragorn.txt'.format(c=cmd,i=infile)
    tmp = subprocess.check_output(cmd, shell=True)
    df = tools.read_aragorn('/tmp/aragorn.txt')
    return df

def default_databases():
    """default blast db table"""

    path = os.path.join(config_path, 'blast_dbs.csv')
    df = {'sprot':{'filename':'sprot.fa','evalue':1e-10},
            'IS':{'filename':'IS.fa','evalue':1e-30},
            'amr':{'filename':'amr.fa','evalue':1e-300}}
    return

def run_annotation(infile, prefix=None, ident=70, threads=4, **kwargs):
    """
    Annotate nucelotide sequences (usually a draft assembly with contigs)
    using prodigal and blast to prokka seqs. Writes a genbank file to the
    same folder.
    Args:
        prefix: prefix for locus_tags
        infile: input fasta file
        outfile: output genbank
        hmmer: run hmmer
    returns:
        a list of SeqRecords with the features
    """

    #get simple name for contig
    def get_contig(x):
        return ('_').join(x.split('_')[:-1])
    if prefix == None:
        prefix = create_locus_tag(infile)
    dbs = ['IS','amr','sprot']
    evalues = [1e-10,1e-100,1e-4]
    #run prodigal
    resfile = prodigal(infile)
    #read in prodigal fasta to dataframe
    df = tools.fasta_to_dataframe(resfile)
    df[['start','end','strand']] = df.description.apply(get_prodigal_coords,1)
    df['feat_type'] = 'CDS'
    df['contig'] = df['name'].apply(get_contig)

    #get target seqs
    seqs = list(SeqIO.parse(resfile,'fasta'))
    #read input file nucleotide seqs
    contigs = SeqIO.to_dict(SeqIO.parse(infile,'fasta'))
    #print (df[:5])
    res = []
    i=0
    for db in dbs:
        fetch_sequence_from_url(db, path=prokkadbdir)
        #make blast db of prokka proteins
        dbname = os.path.join(prokkadbdir,'%s.fa' %db)
        tools.make_blast_database(dbname, dbtype='prot')
        print ('blasting %s ORFs to %s' %(len(seqs),db))
        bl = tools.blast_sequences(dbname, seqs, maxseqs=1, evalue=evalues[i],
                                    cmd='blastp', show_cmd=True, threads=threads, **kwargs)
        bl = bl[bl.pident>ident]
        if len(bl)==0:
            i+=1
            continue
        bl[['protein_id','gene','product','cog']] = bl.stitle.apply(prokka_header_info,1)

        cols = ['qseqid','sseqid','pident','sstart','send','protein_id','gene','product']
        #print (len(bl))
        bl = bl.sort_values(['qseqid','pident'], ascending=False).drop_duplicates(['qseqid'])[cols]
        #print (len(bl))

        #merge blast result with prodigal sequences
        found = df.merge(bl,left_on='name',right_on='qseqid',how='right')
        #get remaining sequences with no hits to this db
        df = df[~df.name.isin(bl.qseqid)]
        #new sequences to blast in the next iteration
        seqs = tools.dataframe_to_seqrecords(df, idkey='name')
        #print (found)
        res.append(found)
        print ('%s sequences unassigned' %len(df))
        i+=1
    #all results together
    res = pd.concat(res)
    #print (res)

    #-------------------------------------------------------
    #run hmmer on unassigned
    print ('running hmmer')
    #write unknowns out
    SeqIO.write(seqs,'unknowns.fa','fasta')
    hmmdf = hmmer('unknowns.fa', threads=threads)
    #print (hmmdf.dtypes)
    res = pd.concat([res,hmmdf], sort=False)

    #get tRNAs with aragorn
    print ('running aragorn')
    arag = aragorn(infile)
    #print (arag)
    res = pd.concat([res,arag], sort=False)

    #remaining unknowns are hypothetical proteins
    unknown = df[~df.name.isin(res.name)]
    unknown['product'] = 'hypothetical protein'
    res = pd.concat([res,unknown], sort=False)

    #post process dataframe
    #res = res.fillna('')
    res['translation'] = res.sequence
    #res['gene'] = res.gene.fillna('')
    res = res.reset_index(drop=True)

    #print (res['product'].value_counts())
    #print (res.dtypes)

    #-------------------------------------------------------
    #we then write all found sequences to seqrecord/features
    l=1  #counter for assigning locus tags
    recs = []
    #group by contig and get features for each protein found
    for c,df in res.groupby('contig'):
        contig = get_contig(c)
        #truncated label for writing to genbank
        label = ('_').join(c.split('_')[:2])
        print (c, len(df), label)
        nucseq = contigs[c].seq
        rec = SeqRecord(nucseq)
        rec.seq.alphabet = generic_dna
        rec.id = label
        rec.name = label
        rec.COMMENT = 'annotated with pygenefinder'
        df = df.sort_values('start')
        qcols = ['gene','product','locus_tag','translation']
        for i,row in df.iterrows():
            row['locus_tag'] = '{p}_{l:05d}'.format(p=prefix,l=l)
            row = row.dropna()
            cols = [c for c in qcols if c in row.index]
            quals = row[cols].to_dict()
            #print (quals)
            feat = SeqFeature(FeatureLocation(row.start,row.end,row.strand), strand=row.strand,
                              type=row.feat_type, qualifiers=quals)
            rec.features.append(feat)
            l+=1
        recs.append(rec)
    return res,recs

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
