"""
    MLST methods for M.bovis.
    Created May 2021
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

import sys, os, string, types, re
import platform, tempfile
import shutil, glob, collections
import itertools
import subprocess
import numpy as np
import pandas as pd
from . import tools, app

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','pathogenie')
module_path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(module_path, 'mlst')

mbovis_scheme = pd.read_csv(os.path.join(datadir,'mbovis_scheme.csv'))
mbovis_db = os.path.join(datadir,'mbovis_db.csv.gz')
ref_proteins = os.path.join(datadir,'Mbovis_AF212297_proteins.fa')

schemes = {'Mbovis-AF212297': mbovis_scheme}

def get_samples_vcf(vcf_file):
    cmd = 'bcftools query -l %s' %vcf_file
    tmp = subprocess.check_output(cmd, shell=True)
    return tmp.decode().split('\n')

def get_nucleotide_sequences(gb_file,out_file,idkey='locus_tag'):
    """protein nucleotide seqs from genbank"""

    recs = SeqIO.to_dict(SeqIO.parse(gb_file,'genbank'))
    chroms = list(recs.keys())
    result = []
    for chrom in chroms:
        rec = recs[chrom]
        for f in rec.features[1:]:
            q=f.qualifiers
            if f.type != 'CDS':
                continue
            seq = rec.seq[f.location.start:f.location.end]
            try:
                new = SeqRecord(seq,id=q[idkey][0])
                result.append(new)
            except:
                #print (q)
                pass
    SeqIO.write(result,out_file,format='fasta')
    return result

def get_consensus(vcf_file, sample, out_file='consensus.fa'):
    """Get consensus sequence from vcf"""

    cmd='bcftools index -f %s' %vcf_file
    subprocess.check_output(cmd, shell=True)
    cmd='cat {r} | bcftools consensus -s {s} {v} > {o}'.format(r=app.mbovis_genome,
                                                v=vcf_file,s=sample,o=out_file)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def diff_profiles(s1, s2):
    return sum(1 for a, b in zip(s1, s2) if a != b)

def find_alleles(fastafile):
    """Find allele by simple matches to the reference table of known sequences.
    Returns:
        dataframe with allele number for each gene
        dataframe with new alleles to add to db
    """

    db = pd.read_csv(mbovis_db)
    names = db.name.unique()
    df = tools.fasta_to_dataframe(fastafile).reset_index()

    result=[]
    new=[]
    for name in names:
        #print (name)
        s = db[db.name==name]
        gene = df[df.name==name]
        #print (gene)
        if len(gene)==0:
            #print (name)
            #missing gene in target
            result.append((name,0))
            continue
        target = gene.iloc[0].sequence
        found = s[s.sequence==target]
        #print (target,found)
        if len(found)>0:
            found = found.iloc[0]
            result.append((name,found.allele))
        else:
            #assign new allele
            newallele = s.allele.max()+1
            result.append((name,newallele))
            new.append([name,newallele,target])
    prof = pd.DataFrame(result,columns=['name','allele'])
    prof['allele'] = prof.allele.astype(int)
    #new additions
    new = pd.DataFrame(new,columns=['name','allele','sequence'])
    return prof, new

def update_mlst_db(new):
    """Update the database of MLST profiles"""

    db = pd.read_csv(mbovis_db)
    db = pd.concat([db,new])
    db.to_csv(mbovis_db, index=False, compression='gzip')
    print ('added %s new alleles' %len(new))
    return

def type_sample(fastafile, outfile, threads=4, overwrite=False):
    """Type a single sample using wgMLST.
    Args:
        fastafile: fasta file to type from assembly or other

        path: output folder for annotations
    Returns:
        dataframe of MLST profile
    """


    if overwrite == True or not os.path.exists(outfile):
        #annotate
        featdf,recs = app.run_annotation(fastafile, threads=threads,
                                        kingdom='bacteria', trusted=ref_proteins)
        #get nucl sequences from annotation
        SeqIO.write(recs,'temp.gb','genbank')
        get_nucleotide_sequences('temp.gb',outfile,idkey='protein_id')

    #find alleles
    res,new = find_alleles(outfile)
    #update db
    update_mlst_db(new)
    return res

def dist_matrix(profiles):
    """Distance matrix of a set of profiles"""

    dist=[]
    for s in profiles:
        x=profiles[s]
        row=[]
        for s in profiles:
            d = diff_profiles(x,profiles[s])
            row.append(d)
        dist.append(row)
    D = pd.DataFrame(dist,columns=profiles.keys(),index=profiles.keys())
    return D

def tree_from_distmatrix(D):
    """tree from distance matrix"""

    from skbio import DistanceMatrix
    from skbio.tree import nj
    ids = list(D.index)
    dm = DistanceMatrix(D.values, ids)
    tree = nj(dm)
    #print(tree.ascii_art())
    return tree

def run_samples(vcf_file, outdir, omit=[], **kwargs):
    """Run samples in a vcf file.
    Args:
        vcf_file: multi sample variant file from previous calling
        outdir: folder for writing intermediate files
    Returns:
        dict of mst profiles
    """

    profs = {}
    samplenames = get_samples_vcf(vcf_file)
    for s in samplenames:
        print (s)
        if s in omit:
            continue
        get_consensus(vcf_file, s)
        outfile = os.path.join(outdir, '%s.fa' %s)
        profile = type_sample('consensus.fa', outfile, **kwargs)
        profs[s] = list(profile.allele)
    return profs

def test():

    return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='pathogenie wgMLST tool. https://github.com/dmnfarrell/pathogenie')
    parser.add_argument("-i", "--input", dest="vcf_file", default=None,
                        help="input vcf from variant calling", metavar="FILE")
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Folder to output intermediate results", metavar="FILE")
    parser.add_argument("-S", "--species", dest="species", default=None,
                        help="set the species")
    parser.add_argument("-t", "--threads", dest="threads", default=None,
                        help="cpu threads to use")
    parser.add_argument("-x", "--test", dest="test",  action="store_true",
                        default=False, help="Test run")

    args = vars(parser.parse_args())
    if args['test'] == True:
        test()
    elif args['vcf_file'] != None:
        profs = run_samples(args['vcf_file'], outdir=args['outdir'])
        D = dist_matrix(profs)
        D.to_csv('dist_mlst.csv',index=False)
        tree = tree_from_distmatrix(D)
        print()
        print(tree.ascii_art())
        tree.write('mlst.newick', 'newick')

if __name__ == '__main__':
    main()
