#!/usr/bin/env python

"""
    pathogenie cmd line tool.
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
import platform
import urllib, hashlib, shutil
import tempfile
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from Bio.Alphabet import generic_dna
from . import tools

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','pathogenie')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
dbdir = os.path.join(config_path, 'db')
customdbdir = os.path.join(config_path, 'custom')
hmmdir = os.path.join(config_path, 'hmms')
prokkadbdir = os.path.join(config_path, 'prokka')
trustedproteindir = os.path.join(config_path, 'proteins')
tempdir = os.path.join(config_path, 'tmp')

db_names = ['card','resfinder','argannot','ncbi','plasmidfinder','ecoh','vfdb',
            'bacteria.16SrRNA','bacteria.23SrRNA']
prokka_db_names = ['sprot_bacteria','sprot_viruses','IS','AMR']
links = {'card':'https://github.com/tseemann/abricate/raw/master/db/card/sequences',
        'resfinder':'https://raw.githubusercontent.com/tseemann/abricate/master/db/resfinder/sequences',
        'vfdb':'https://raw.githubusercontent.com/tseemann/abricate/master/db/vfdb/sequences',
        'ncbi':'https://raw.githubusercontent.com/tseemann/abricate/master/db/ncbi/sequences',
        'argannot':'https://raw.githubusercontent.com/tseemann/abricate/master/db/argannot/sequences',
        'ecoh':'https://raw.githubusercontent.com/tseemann/abricate/master/db/ecoh/sequences',
        'plasmidfinder':'https://raw.githubusercontent.com/tseemann/abricate/master/db/plasmidfinder/sequences',
        'sprot_bacteria':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/sprot',
        'sprot_viruses':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Viruses/sprot',
        'sprot_archaea':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Archaea/sprot',
        'amr':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/AMR',
        'IS':'https://raw.githubusercontent.com/tseemann/prokka/master/db/kingdom/Bacteria/IS',
        'bacteria.16SrRNA': 'https://raw.githubusercontent.com/dmnfarrell/pathogenie/master/db/bacteria.16SrRNA.fna',
        'bacteria.23SrRNA': 'https://raw.githubusercontent.com/dmnfarrell/pathogenie/master/db/bacteria.23SrRNA.fna',
        'HAMAP.hmm':'https://github.com/tseemann/prokka/raw/master/db/hmm/HAMAP.hmm',
        'TIGRFAMS':'ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz',
        'PFAM-A':'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz' }

hmm = 'https://github.com/tseemann/prokka/raw/master/db/hmm/HAMAP.hmm'

if not os.path.exists(config_path):
    try:
        os.makedirs(config_path, exist_ok=True)
    except:
        os.makedirs(config_path)

os.makedirs(tempdir, exist_ok=True)

def check_platform():
    """See if we are running in Windows"""

    if platform.system() == 'Windows':
        print('checking binaries are present')
        tools.fetch_binaries()
    return

check_platform()

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
    """Create a genbank style locus tag"""

    name = os.path.basename(filename)
    x = hashlib.md5(name.encode('utf-8')).hexdigest()
    x = ''.join(i for i in x if not i.isdigit())
    x = x.upper()[:7]
    return x

def prodigal(infile):
    """Run prodigal"""

    cmd = tools.get_cmd('prodigal')
    #name = os.path.splitext(infile)[0]
    name = os.path.join(tempdir,'prodigal')
    cmd = '{c} -i {i} -c -m -a {n}.faa -f gff -o {n}.gff -p single'.format(i=infile,c=cmd,n=name)
    subprocess.check_output(cmd, shell=True)
    resfile = os.path.join(tempdir, 'prodigal.faa')
    return resfile

def get_prodigal_coords(x):
    """Extract prodigal coordinate information"""

    s = re.split('\#|\s',x.replace(' ',''))
    coords = [int(i) for i in s[1:4]]
    coords[0] = coords[0]-1
    return  pd.Series(coords)

def prokka_header_info(x):
    s = re.split('~~~',x)
    return pd.Series(s)

def uniprot_header_info(x):
    """Extract uniprot header info to series"""

    s = x.split('|')
    y = re.split('OS=|GN=|OX=|PE=|SV=',s[2])
    product = y[0]
    pid = s[1]
    gene = y[3]
    return pd.Series([pid,gene,product,''])

def hmmer(infile, threads=4, hmm_file=None):
    """Run hmmer"""

    def get_contig(x):
        return ('_').join(x.split('_')[:-1])

    hmmpresscmd = tools.get_cmd('hmmpress')
    hmmscancmd = tools.get_cmd('hmmscan')

    df = tools.fasta_to_dataframe(infile)
    fetch_sequence_from_url('HAMAP.hmm', hmmdir, ext='')
    db = os.path.join(hmmdir,'HAMAP.hmm')
    #press the database
    cmd = '%s -f %s' %(hmmpresscmd,db)
    tmp = subprocess.check_output(cmd, shell=True)
    outfile = os.path.join(tempdir, 'hmm.txt')
    cmd = "{c} --noali --notextw --acc -E 1e-4 --cpu {t} --tblout {o} -o hmm.out {db} {i}".format(
                                c=hmmscancmd,t=threads,db=db,i=infile,o=outfile)
    print (cmd)
    tmp = subprocess.check_output(cmd, shell=True)
    h = tools.read_hmmer3(outfile)
    if h is None:
        return
    df = h.merge(df,on='name',how='left')
    #print (df[:5])
    #get coords from prodigal fasta heading if available
    df[['start','end','strand']] = df.description.apply(get_prodigal_coords,1)
    df['contig'] = df['name'].apply(get_contig)
    #remove temp file
    os.remove(outfile)
    return df

def aragorn(infile):
    """Run aragorn"""

    outfile = os.path.join(tempdir, 'aragorn.txt')
    cmd = tools.get_cmd('aragorn')
    cmd = '{c} -l -gcbact -t -w {i} -o {o}'.format(c=cmd,i=infile,o=outfile)
    tmp = subprocess.check_output(cmd, shell=True)
    df = tools.read_aragorn(outfile)
    return df

def default_databases():
    """default blast db table"""

    path = os.path.join(config_path, 'blast_dbs.csv')
    df = {'sprot_bacteria':{'filename':'sprot_bacteria.fa','evalue':1e-10},
          'sprot_viruses':{'filename':'sprot_viruses.fa','evalue':1e-10},
            'IS':{'filename':'IS.fa','evalue':1e-30},
            'amr':{'filename':'amr.fa','evalue':1e-300}}
    return

def get_files_in_path(path):
    """Find names of custom sequence files"""

    return glob.glob(path+'/*')

def add_protein_db(filename):
    """Add preferred protein sequences to be used for annotation"""

    path = trustedproteindir
    os.makedirs(path, exist_ok=True)
    newpath = os.path.join(path, os.path.basename(filename))
    shutil.copy(filename, newpath)
    return

def run_annotation(infile, prefix=None, ident=70, threads=4,
                   kingdom='bacteria', trusted=None, trusted_format='uniprot',
                   callback=None,
                   **kwargs):
    """
    Annotate nucelotide sequences (usually a draft assembly with contigs)
    using prodigal and blast to prokka seqs. Writes a genbank file to the
    same folder.
    Args:
        infile: input fasta file
        prefix: optional prefix for locus_tags
        ident: minimum percent identity for blast results, default 70
        threads: cpu threads to use
        kingdom: 'bacteria', 'viruses' or 'archaea'
        trusted: optional fasta file of trusted protein sequences for blast
    returns:
        a dataframe of annotations and a list of SeqRecords with the features, one per contig
    """

    #get simple name for contig
    def get_contig(x):
        return ('_').join(x.split('_')[:-1])
    if prefix == None:
        prefix = create_locus_tag(infile)
    sprot = 'sprot_%s' %kingdom
    dbs = ['IS','amr',trusted,sprot]
    evalues = [1e-10,1e-100,1e-4,1e-4]
    #run prodigal
    resfile = prodigal(infile)
    #read in prodigal fasta to dataframe
    df = tools.fasta_to_dataframe(resfile)
    df[['start','end','strand']] = df.description.apply(get_prodigal_coords,1)
    df['feat_type'] = 'CDS'
    df['contig'] = df['name'].apply(get_contig)
    df['sequence'] = df.sequence.str.rstrip('*')
    #get target seqs
    seqs = list(SeqIO.parse(resfile,'fasta'))
    #remove temp prodigal file
    os.remove(resfile)
    #remove trailing asterisks
    #seqs = [s.rstrip("*") for s in seqs]
    #read input file nucleotide seqs
    contigs = SeqIO.to_dict(SeqIO.parse(infile,'fasta'))
    #print (df[:5])
    res = []
    i=0
    for db in dbs:
        if db is None:
            i+=1
            continue
        fetch_sequence_from_url(db, path=prokkadbdir)
        #make blast db of prokka proteins
        if db in ['IS','amr',sprot]:
            dbname = os.path.join(prokkadbdir,'%s.fa' %db)
        else:
            dbname = db
        tools.make_blast_database(dbname, dbtype='prot')
        print ('blasting %s ORFs to %s' %(len(seqs),db))
        bl = tools.blast_sequences(dbname, seqs, maxseqs=1, evalue=evalues[i],
                                    cmd='blastp', show_cmd=True, threads=threads, **kwargs)
        bl = bl[bl.pident>ident]
        #print (bl)
        if len(bl)==0:
            i+=1
            continue
        #get hit info from header (prokka format)
        try:
            bl[['protein_id','gene','product','cog']] = bl.stitle.apply(prokka_header_info,1)
        except:
            bl[['protein_id','gene','product','cog']] = bl.stitle.apply(uniprot_header_info,1)
        cols = ['qseqid','sseqid','pident','sstart','send','protein_id','gene','product']

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
    if kingdom == 'bacteria':
        print ('running hmmer')
        #write unknowns out
        SeqIO.write(seqs,'unknowns.fa','fasta')
        hmmdf = hmmer('unknowns.fa', threads=threads)
        if hmmdf is not None:
            res = pd.concat([res,hmmdf], sort=False)

    #get tRNAs with aragorn
    print ('running aragorn')
    arag = aragorn(infile)
    #print (arag)
    res = pd.concat([res,arag], sort=False)

    #remaining unknowns are hypothetical proteins
    unknown = df[~df.name.isin(res.name)].copy()
    unknown['product'] = 'hypothetical protein'
    res = pd.concat([res,unknown], sort=False)

    #post process dataframe
    #res = res.fillna('')
    res['translation'] = res.sequence
    res['length'] = res.sequence.str.len()
    #res['gene'] = res.gene.fillna('')
    res = res.reset_index(drop=True)

    #print (res['product'].value_counts())
    #print (res.dtypes)

    #-------------------------------------------------------
    #we then write all found sequences to seqrecord/features
    l=1  #counter for assigning locus tags
    recs = []
    #print (res.iloc[])
    #group by contig and get features for each protein found
    for c,df in res.groupby('contig'):
        contig = get_contig(c)
        #truncated label for writing to genbank
        label = ('_').join(c.split('_')[:2])
        #print (c, len(df), label)
        nucseq = contigs[c].seq
        rec = SeqRecord(nucseq, annotations={"molecule_type": "DNA"})
        #rec.seq.alphabet = generic_dna
        rec.id = label
        rec.name = label
        rec.COMMENT = 'annotated with pathogenie'
        df = df.sort_values('start')
        qcols = ['gene','product','locus_tag','translation','length']
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
    print ('done')
    return res,recs

def annotate_files(recs, keys=None, outdir='annot', kingdom='bacteria',
                    overwrite=False):
    """Annotate a set of nucleotide seqrecords.
    Args:
        recs: list of SeqRecords
        outdir: output folder to save genbank files
        kingdom: 'bacteria', 'viruses' or 'archaea'
    Returns: a dataframe with annotations.
    """

    res = []
    for label in recs:
        if keys!=None and label not in keys:
            continue
        rec = recs[label]
        gbfile = os.path.join(outdir,label+'.gbk')
        if os.path.exists(gbfile) and overwrite == False:
            featdf = tools.genbank_to_dataframe(gbfile)
            featdf['sequence'] = featdf.translation
        else:
            seq = recs[label]
            filename = os.path.join('temp',label+'.fasta')
            SeqIO.write(seq,filename,'fasta')
            featdf,protrecs = run_annotation(filename, threads=10, kingdom=kingdom)
            tools.recs_to_genbank(protrecs, gbfile)
        featdf['label'] = label
        res.append(featdf)
    res = pd.concat(res)
    return res

def get_similar_sequences(protname, annot):
    """Extract similar sequences by name from a set of annotations"""

    seqs = []
    for i,df in annot.groupby('label'):
        s = df[df['product']==protname]
        if len(s)==0:
            continue
        s = s.iloc[0]
        seq = SeqRecord(Seq(s.sequence),id=s.label)#,description=s.host)
        seqs.append(seq)
    return seqs

def find_mutations(recs, ref):
    """Find the mutations in a set of protein records relative to a
    reference protein sequence"""

    mutations = {}
    positions = []
    for rec in recs:
        aln = tools.clustal_alignment(seqs=[ref, rec])
        #print (aln)
        x = []
        for pos in range(len(aln[0])):
            refaa = aln[0,pos]
            aa = aln[1,pos]
            if aa != refaa and aa!='-':
                #print (refaa, aln[:,pos], aa)
                mut = refaa+str(pos+1)+aa
                x.append(mut)
        if len(x)>0:
            mutations[rec.seq] = x
    return mutations

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
    parser = ArgumentParser(description='pathogenie command line')

    parser.add_argument("-f", "--fasta", dest="filenames", nargs='*',
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-p", "--path", dest="path",
                        help="input fasta file", metavar="FILE")
    parser.add_argument("-o", "--out", dest="outdir", default='results_pygf',
                        help="output folder", metavar="FILE")
    parser.add_argument("-k", "--kingdom", dest="kingdom", default='bacteria',
                        help="kingdom for annotation")
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
