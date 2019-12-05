"""
python module with various methods for bacterial annotation and comparative analysis
"""

from __future__ import print_function
import sys,os,subprocess,glob,shutil
from Bio import Entrez
Entrez.email = 'A.N.Other@example.com'
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import matplotlib as mpl
import pylab as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, LogNorm
import numpy as np
import pandas as pd

def fastq_to_dataframe(f, size=None):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    import HTSeq
    ext = os.path.splitext(f)[1]
    if ext=='.fastq' or ext=='.gz':
        ffile = HTSeq.FastqReader(f, "solexa")
    elif ext == '.fa':
        ffile = HTSeq.FastaReader(f)
    else:
        return
    if size != None:
        sequences = [(s.name, s.seq, s.descr) for s in islice(fastfile, i, i+size)]
    else:
        sequences = [(s.name,s.seq) for s in ffile]
    df = pd.DataFrame(sequences,columns=['id','seq'])
    return df

def read_length_dist(df):

    df['length'] = df.seq.str.len()
    bins = np.linspace(1,df.length.max(),df.length.max())
    x = np.histogram(df.length,bins=bins)
    return x

def dataframe_to_fasta(df, seqkey='translation', idkey='locus_tag',
                     descrkey='description',
                     outfile='out.faa'):
    """Genbank features to fasta file"""

    seqs=[]
    for i,row in df.iterrows():
        if descrkey in df.columns:
            d=row[descrkey]
        else:
            d=''
        rec = SeqRecord(Seq(row[seqkey]),id=row[idkey],
                            description=d)
        seqs.append(rec)
    SeqIO.write(seqs, outfile, "fasta")
    return outfile

def fasta_to_dataframe(infile, header_sep=None, key='name', seqkey='sequence'):
    """Get fasta proteins into dataframe"""

    recs = SeqIO.parse(infile,'fasta')
    keys = [key,seqkey,'description']
    data = [(r.name,str(r.seq),str(r.description)) for r in recs]
    df = pd.DataFrame(data,columns=(keys))
    df['type'] = 'CDS'
    #fix bad names
    if header_sep not in ['',None]:
        df[key] = df[key].apply(lambda x: x.split(header_sep)[0],1)
    df[key] = df[key].str.replace('|','_')
    return df

def run_blastn(database, query, params='-e .1 -G 10'):
    """Run blast"""

    out = 'blast_result.csv'
    cmd = 'blastall -d %s -i %s -p blastn -m 8 %s > %s' %(database,query,params,out)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    #zipfile(out+'.xml', remove=True)
    return

def local_blast(fastafile, database, ident=100, params='-e 10 -a 2', results='all'):
    """Blast a blastdb and save hits to a dataframe.
       Args:
            fastafile: file with queries
            database: name of blast db
            ident: cutoff percentage identity
            params: custom blast parameters (see blastall -h)
            results: return 'all' or 'best'
        Returns: dataframe of hits"""

    outname = 'blast_result.csv'
    run_blastn(database, fastafile, params)
    cols = ['query', 'subj', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore']
    res = pd.read_csv(outname, names=cols, sep='\t')
    print ('found %s hits in db' %len(res))
    res = res[res['pident']>=ident]
    if results == 'best':
        res = res.groupby('query').first().reset_index()
    queries = fasta_to_dataframe(fastafile).reset_index()
    res = res.merge(queries, left_on='query', right_on='name', how='left')
    print ()
    return res

def clustal_alignment(filename=None, seqs=None, command="clustalw"):
    """Align 2 sequences with clustal"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline(command, infile=filename)
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', 'clustal')
    return align

def muscle_alignment(filename=None, seqs=None):
    """Align 2 sequences with muscle"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=name+'.txt')
    stdout, stderr = cline()
    align = AlignIO.read(name+'.txt', 'fasta')
    return align

def align_nucmer(file1, file2):
    cmd='nucmer --maxgap=500 --mincluster=100 --coords -p nucmer %s %s' %(file1, file2)
    print (cmd)
    subprocess.check_output(cmd,shell=True)
    df = read_nucmer_coords('nucmer.coords')
    return df

def read_nucmer_coords(cfile):
    cols=['S1','E1','S2','E2','LEN 1','LEN 2','IDENT','TAG1','TAG2']
    a=pd.read_csv(cfile,sep='[\s|]+',skiprows=5,names=cols,engine='python')
    a = a.sort_values(by='TAG2',ascending=False)
    return a

def align_reads(file1, file2, idx, out):
	"""align reads to ref"""

    #cmd = 'bowtie2 -x %s -1 %s -2 %s --threads 6 | samtools view -bS - > %s' %(idx,files[0],files[1],out)
	cmd = 'bwa mem -M -t 8 %s %s %s | samtools view -bS - > %s' %(idx,files[0],files[1],out)
	if not os.path.exists(out):
		print (cmd )
		subprocess.check_output(cmd, shell=True)
	return

def align_info(bamfile):
    cmd = 'samtools flagstat %s' %bamfile
    temp=subprocess.check_output(cmd, shell=True)
    print (temp)
    return

def variants_call(name, ref, out):
    bamfile = '%s/%s.bam' %(out,name)
    cmd = 'samtools sort {b} > {b}.sorted && samtools index {b}.sorted'.format(b=bamfile)
    print (cmd)
    #subprocess.check_output(cmd, shell=True)
    cmd = 'samtools mpileup -uf genomes/{r}.fa {b}.sorted | bcftools call -mv \
    > {o}/{n}.vcf'.format(b=bamfile,n=name,r=ref,o=out)
    print (cmd)
    #subprocess.check_output(cmd, shell=True)
    cmd = 'bedtools intersect -a {gff} -b {o}/{n}.vcf -wa -u > {o}/{n}_variants.bed'.format(n=name,r=ref,gff=gff,o=out)
    print (cmd)

def search_genbank(term='', filt=None):
    request = Entrez.esearch(db="nuccore", term=term, field="title", FILT=filt, rettype='xml')
    result = Entrez.read(request)
    idlist = result['IdList']
    return idlist

def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("nucleotide",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print ("An error occurred while retrieving the annotations.")
        print ("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)
    print ("Retrieved %d annotations for %d genes" % (len(annotations),
            len(id_list)))
    return annotations

def retrieve_sequences(id_list):
    """get entrez sequence"""

    request = Entrez.epost("nucleotide",id=",".join(id_list))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="nucleotide",retmode="xml", webenv=webEnv, query_key=queryKey)
    recs={}
    for r in Entrez.parse(handle):
        recs[r['GBSeq_primary-accession']] = r
    return recs

def recs_to_fasta(recs,outfile):
    res=[]
    for i in recs:
        r=recs[i]
        res.append( [r['GBSeq_primary-accession'],r['GBSeq_sequence'],r['GBSeq_definition']] )
    df=pd.DataFrame(res,columns=['id','sequence','description'])
    dataframe_to_fasta(df,outfile=outfile,seqkey='sequence',idkey='id')
    return

def get_gilist(accs):
    query  = " ".join(accs)
    handle = Entrez.esearch(db="nucleotide",term=query,retmax=10000)
    gilist = Entrez.read(handle)['IdList']
    return gilist

def read_blast_nr(filename):
    blastcols = ['contig','gid','name','accession','pident','length',
                 '?', '61', 'qstart', 'qend', 'sstart', 'send',
                  'evalue', 'score']
    bl = pd.read_csv('909D3A_blast_nr.csv',names=blastcols)
    #print bl[bl.contig=='NODE_1_length_485821_cov_65.8942']
    g = bl.groupby(['contig','accession']).agg({'score':np.sum,'length':np.sum,'pident':np.max})
    g = g.sort_values('score',ascending=False).reset_index()
    return g

#prokka and roary utils

def get_gene_name(x):
    p=str(x['product'])
    if x.gene is np.nan:
        s=p.split(':')[0][:35]
        if s == None:
            s=p.split(';')[0][:35]
        return s
    return x.gene

def get_aro(x):
    x=str(x)
    if 'ARO' in x:
        s=x[x.find("[")+1:x.find("]")].split(';')[0]
        return s

def get_product(x):
    x=str(x)
    s = x.split(';')[0]
    #print s
    return s

def prokka_results(path,names):
    """parse prokka results for multiple files"""
    res=[]
    for n in names:
        f = '%s/%s/%s.tsv' %(path,n,n)
        df=pd.read_csv(f,sep='\t')
        df['isolate'] = n
        print (n, len(df))
        #print df
        res.append(df)
    res=pd.concat(res)
    res['gene'] = res.apply( lambda x: get_gene_name(x),1)
    res['cat'] = res['product'].apply(lambda x: apply_cat(x))
    #res['aro'] = res['product'].apply(lambda x: get_aro(x))
    res['fam'] = res['product'].apply(lambda x: get_product(x))
    return res

def plot_products(res):
    #res=res.merge(info,left_on='isolate',right_on='id')
    g = res.groupby('isolate').agg({'product':np.size})#.reset_index()
    g.plot(kind='bar',figsize=(10,5),legend=False)
    plt.title('Annotated Proteins',fontsize=20)
    plt.ylabel('proteins')
    plt.tight_layout()
    plt.savefig('product_counts_ecoli.png')
    return

def get_presence_absence(df, cols=None):
    """parse roary file"""
    if cols is None:
        cols = df.columns[14:]
    x=df.copy()
    x['cat'] = x.Annotation.apply(lambda x: apply_cat(x))
    x = x.set_index(['cat','Gene','Annotation'])
    x = x[cols].notnull().astype('int')
    #x = x.loc[x.index.dropna()]
    return x

def apply_cat(x):
    keys=['ARO','efflux','adhesin','LEE','porin','stress',
          'secretion system', 'bacteriophage',
          'membrane','prophage','secreted','IS','insertion','transposase','integrase',
          'virulence','protease','stress','toxic','phage','kinase','phosphatase',
          'hypothetical','membrane','binding','rna','ribosomal','tRNA','methyltransferase',
         'polymerase','DNA','Transcription','lipoprotein','protease']
    for i in keys:
        if x is np.nan: return
        if i in x:
            return i
    return 'other'

def genes_clustermap(x,xticklabels=0,title=''):
    """plot cluster map of genes"""
    #x = x[x.sum(1)>=1]
    sys.setrecursionlimit(20000)
    clrs = ["lightgray", "blue"]
    t_cmap = ListedColormap(sns.color_palette(clrs).as_hex())
    if len(x)>50:
        yticklabels=0
    if len(x.T)>50:
        xticklabels=0
    cg=sns.clustermap(x,xticklabels=xticklabels,yticklabels=1,cmap=t_cmap,figsize=(12,7))
    cg.fig.suptitle(title)
    cg.fig.subplots_adjust(right=0.8)
    return

def get_roary_gene(gene, path='compare_acinet/roary'):
    x = pd.read_csv(os.path.join(path, 'gene_presence_absence.csv'))

    if not gene in list(x.Gene):
        print ('no such gene name')
        return
    r = x[x.Gene==gene].iloc[0].dropna()

    samples = x.columns[14:]
    lbls = r.loc[samples].to_dict()
    lbls = dict((lbls[k], k) for k in lbls)
    f = '%s/pan_genome_sequences/%s.fa.aln' %(path,gene)
    if not os.path.exists(f):
        return
    aln = AlignIO.read(f,'fasta')
    for a in aln:
        a.id = lbls[a.id]
    return aln

def get_prot_seq(name,gene):
    #get prot seqs from annot fasta files
    f='annot_scaff/{n}/{n}.faa'.format(n=name)
    #print f
    seqs = SeqIO.to_dict(SeqIO.parse(f,'fasta'))
    return seqs[gene]

def get_roary_protein(gene, path='roary', samples=None):
    #get protein seqs for a gene in all samples
    x = pd.read_csv(os.path.join(path, 'gene_presence_absence.csv'))
    if samples is None:
        samples = x.columns[14:]
    #print x[x.Gene==gene]
    r = x[x.Gene==gene].iloc[0].dropna()
    lbls = r.loc[samples].to_dict()
    lbls = dict((lbls[k], k) for k in lbls)
    seqs=[]
    for a in samples:
        try:
            s = get_prot_seq(a,r[a])
            #print s
        except:
            continue
        seqs.append(s)
    return seqs, lbls

def draw_features(rec):
    from dna_features_viewer import BiopythonTranslator
    graphic_record = BiopythonTranslator().translate_record(rec)
    ax, _ = graphic_record.plot(figure_width=20)
    plt.title(rec.id)
    plt.show()

css='''
body {
    font-family: arial;
}
h1,h2,h3,h4,h5,h6 {margin:0.2em 0 0.25em 0; display:block;
  font-family: Arial}

@font-face{
  font-family: open sans;
  src: url('OpenSans-Regular.ttf');
}

.tinytable {
    font-family:monospace;
    font-size:12px;
}
.tinytable table{
    margin:2px;
    padding:3;
}.tinytable td{
    border-spacing: 0;
    border:0px solid #b2a7a7;
    border-width:0px 0px 0px 0px;
    width: 80px;
    white-space: nowrap;
}
.tinytable th{
    text-align: left;
    background-color: #F1F1F1;
    border-width:0px 0px 0px 0px;
}
'''

def contigs_from_prokka(gff_file, key):
    """get contig features from prokka gff"""

    from BCBio import GFF
    res=[]
    for rec in GFF.parse(open(gff_file)):
        found=False
        x=[]
        if len(rec.seq)>80000:
            continue
        for feat in rec.features:
            q=feat.qualifiers
            #print q
            try:
                p = q['product'][0]
                q['source'] = p
                tag=q['locus_tag'][0]
            except:
                continue
            if key in p or key=='any':
                found=True
            x.append([rec.id,p,str(feat.location),tag,len(rec.seq)])
        if found is True:
            res.extend(x)
            #draw_features(rec)
    df= pd.DataFrame(res,columns=['contig','product','loc','protein_id','size'])
    return df

def contig_report(files, outfile, key=None):
    f=open(outfile,'w')
    s='<html><head> <style> %s </style></head>' %css
    s+='<h3>contigs with prokka annotated features</<h3>'
    found=[]
    for infile in files:
        df = contigs_from_prokka(infile, key=key)
        n=os.path.basename(infile)
        df['name'] = n
        c = list(df.contig.unique())
        df=df[-df.contig.isin(found)]
        s+='<div>\n'
        if len(df)>0:
            s+='<h4>%s</h4>' %(n)
            s += df.set_index(['size','contig','protein_id']).sort_index().to_html(classes='tinytable')
        s+='</div>\n'
        found.extend(c)
    f.write(s)
    return

#phylo trees

def get_tree(aln, kind='nj'):
    from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return dm, tree

def clear_clades(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if 'Inner' in clade.name :
            clade.name = ''

        names[clade.name] = clade
    return names

def draw_tree(tree, root=None, labels=None, clear=True, title=''):
    f,ax=plt.subplots(figsize=(10,8))
    if clear == True:
        try:
            clear_clades(tree)
        except:
            pass
    if root != None:
        tree.root_with_outgroup(root)
    if labels != None:
        for clade in tree.get_terminals():
            key = clade.name
            if key in labels:
                clade.name = '%s; %s' %(key,labels[key])
                #clade.name = labels[key]

    Phylo.draw(tree,axes=ax,axis=('off',),branch_labels=None,show_confidence=False)
    #Phylo.draw_graphviz(tree,axes=ax,node_size=0)
    ax.set_title(title,fontsize=16)
    return f, tree

def ml_tree(aln, name):
    from Bio.Phylo.Applications import PhymlCommandline
    AlignIO.write(aln, '%s.phy' %name, 'phylip-relaxed')
    cmdline = PhymlCommandline(input='%s.phy' %name, datatype='nt', alpha='e', bootstrap=100)
    print (cmdline)
    cmdline()
    mtree = Phylo.read('%s.phy_phyml_tree.txt' %name,'newick')
    return mtree

def show_alignment(aln, diff=False, offset=0):
    """
    Show a sequence alignment
        Args:
            aln: alignment
            diff: whether to show differences
    """

    ref = aln[0]
    l = len(aln[0])
    n=60
    chunks = [(i,i+n) for i in range(0, l, n)]
    for c in chunks:
        start,end = c
        lbls = np.arange(start,end,10)-offset
        print (('%-21s' %'name'),''.join([('%-10s' %i) for i in lbls]))
        print (('%21s' %ref.id[:20]), ref.seq[start:end])

        if diff == True:
            for a in aln[1:]:
                diff=''
                for i,j in zip(ref,a):
                    if i != j:
                        diff+=j
                    else:
                        diff+='-'
                name = a.id[:20]
                print (('%21s' %name), diff[start:end])
        else:
            for a in aln[1:]:
                name = a.id[:20]
                print (('%21s' %name), a.seq[start:end])
    return

#ete tree routines

def ete_draw(t,fname=None,title='',mode='r',outfile='ete_tree.png'):
    from ete3 import Tree, TreeStyle, Phyloxml, TextFace
    t.children
    ts = TreeStyle()
    #ts.branch_vertical_margin = 10
    ts.show_scale=False
    ts.scale =  800
    ts.show_leaf_name = False    
    ts.mode = mode
    ts.title.add_face(TextFace(title, fsize=14), column=0)        
    t.render(outfile,dpi=300,h=800,tree_style=ts)
    return ts

def ete_colors(tree, colors=None):
    from ete3 import TreeStyle, TextFace, NodeStyle
    for node in tree.traverse():  
        if node.name in colors:
            #lbl = labels[node.name]
            clr = colors[node.name]            
            #node.add_face(TextFace(node.name, fsize=14), column=0, position='branch-right')
            #f2=TextFace(lbl)
            #node.add_face(f2, column=1, position='aligned')
            #f2.margin_left = 10
            nstyle = NodeStyle()
            nstyle["fgcolor"] = clr
            nstyle["size"] = 15
            node.set_style(nstyle)

def ete_labels(tree, labels, column=1, position='aligned'):
    from ete3 import TreeStyle, TextFace
    for node in tree.traverse():  
        if node.name in labels:
            lbl = labels[node.name]            
            f2=TextFace(lbl)
            node.add_face(f2, column=column, position=position)
            f2.margin_left = 10
            
def get_ete_tree(filename):    
    from ete3 import Tree, Phyloxml
    p = Phyloxml()
    p.build_from_file(filename)   
    t = p.get_phylogeny()[0]
    return t
