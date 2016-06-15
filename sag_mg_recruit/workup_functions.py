import pandas as pd
import gzip
import os.path as op
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
%matplotlib inline

import glob
import os

from pysam import FastxFile


def readfx(fastx):
    fx = ""
    try:
        fx = FastxFile(fastx)
        for f in fx:
            yield f.name, f.sequence, f.quality
    finally:
        if fx:
            fx.close()


gzopen = lambda x: gzip.open(x) if x.endswith(".gz") else open(x)


def plot_read_size(fastq):
    readsizes = defaultdict(lambda: 0)
    readcount = 0
    for name, seq, qual in readfx(fastq):
        readsizes[len(seq)] += 1
        readcount += 1
    hist = pd.DataFrame.from_dict(readsizes, orient='index')
    hist['length']=hist.index
    hist['read_count']=hist[0]
    fig = plt.plot(hist['length'], hist['read_count'], color='b')
    plt.ylabel('# reads')
    plt.xlabel('read length')
    filename = op.dirname(fastq)+"_".join(op.basename(fastq).split(".")[:-1])+".png"
    name = "_".join(op.basename(fastq).split(".")[:-1])+".png"
    plt.title('%s: %s total reads' % (name, readcount))
    return plt, readcount


def get_recruit_info(gcov):
    countfile = gcov.replace("genomecoverage", "aln_count")
    with open(countfile) as infile:
        recruit_count = infile.read().split()[1].strip()
        
    metagenome = op.basename(gcov).split("_vs_")[0].split("_")[0]
    coverage = pd.read_csv(gcov, sep="\t", header=None)
    mean_per_contig = coverage.groupby([0])[2].mean() #.to_dict()
    sum_per_contig = coverage.groupby([0])[2].sum() #.to_dict()
    contig_size = coverage.groupby([0])[1].max()+1
    mean_sag_coverage = mean_per_contig.mean()
    totalbp = contig_size.sum()
    uncovered_bp = sum(coverage[2]==0)
    pct_covered = (totalbp - uncovered_bp)/totalbp * 100
    total_scaffold = len(sum_per_contig)
    uncovered_contig = sum(sum_per_contig==0)
    pct_scaffolds_covered = (total_scaffold - uncovered_contig)/total_scaffold *100
    sag = "_".join(op.basename(gcov).split("_vs_")[1].strip(".genomecoverage").split("_")[:-1])
    cols = ['method','sag','metagenome','Percent_scaffolds_with_any_coverage', 'Percent_of_reference_bases_covered', 'Average_coverage', 'unambiguous']
    data = ['bwa', sag, metagenome, 
           pct_scaffolds_covered,
           pct_covered, 
           mean_sag_coverage,
           recruit_count]
    df = pd.DataFrame(data, index=cols).transpose()
    return df

def genome_cov_table(gcov_list):
    cols = ['sag','metagenome','Percent_scaffolds_with_any_coverage', 'Percent_of_reference_bases_covered', 'Average_coverage', 'unambiguous']
    big = pd.DataFrame(columns=cols)
    for g in gcov_list:
        new=get_recruit_info(g)
        big = pd.concat([big, new])
    return big



    


