from __future__ import print_function
import subprocess
import os
import pandas as pd
from pandas import Series
from itertools import groupby
import click
import glob
import re
import os.path as op
import logging
from sarge import get_stderr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scgc.utils import *
from scgc.fastx import readfx

'''
read recruitment of metagenomes to sags using bbmap
'''


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """read recruitment."""
    pass


def run_flash(prefix, fastq1, fastq2=None, mismatch_density=0.05, min_overlap=35,
        max_overlap=150, cores=1):
    """Joins interleaved FASTQ file using `flash`.

    Args:
        prefix (str): prefix of output files
        fastq (str): file path to fastq
        out_file (str): path of desired gzip compressed fastq
        mismatch_density (Optional[float]): mismatch density of overlapping region
        min_overlap (Optional[int]): minimum expected read overlap
        max_overlap (Optional[int]): maximum expected read overlap
        cores (Optional[int]): threads allocated to flash

    Returns:
        str
    """

    outsuffix = [".extendedFrags.fastq.gz", "out.notCombined_1.fastq.gz", ".notCombined_2.fastq.gz", ".hist", ".histogram"]
    outfiles = [prefix+i for i in outsuffix]
    print("FASTQ1 for FLASH is:", fastq1)
    #with file_transaction(outfiles) as tx:
    with file_transaction(outfiles) as tx_outfiles:
        if fastq2 is None:
            print("only one fastq file provided, assuming it is interleaved", file=sys.stderr)
            cmd = ("flash --interleaved-input -x {mismatch} -m {min} "
                   "-M {max} -t {threads} -z -o {prefix} {fastq1} ").format(mismatch=mismatch_density,
                        min=min_overlap, max=max_overlap, threads=cores, prefix=prefix,
                        fastq1=fastq1)
            print("THE COMMAND WAS:", cmd)
            run(cmd, description="Joining reads with flash")
        else:
            cmd = ("flash -m {min} -M {max} -x {mismatch} -t {threads} "
           "-z -o {prefix} {fastq1} {fastq2}".format(mismatch=mismatch_density,
                        min=min_overlap, max=max_overlap, threads=cores,
                        fastq1=fastq1, fastq2=fastq2, prefix=prefix))
            run(cmd, description="Joining reads with flash")
    return outfiles


def join_stats(inhist, fastq1, fastq2=None, prefix = "", outdir=""):
    '''
    Print join stats and png of read size distribution
    Args:
        inhist (str): path to .hist output created by flash
        prefix (str): location to write output to, defaults to current directory 
    Output:
        joined read size distribution plot as png
        joined read stats as a tab-separated text file
    '''
    name = os.path.basename(inhist).replace(".hist","")

    outname = os.path.join(outdir, "{prefix}_joinstats.txt".format(**locals()))
    png_name = os.path.join(outdir, "{prefix}_joinstats.png".format(**locals()))
    
    hist1 = pd.read_table(inhist, header=None)
    hist1.columns = ["length", "read_count"]
    fig = plt.plot(hist1['length'], hist1['read_count'], color='b')
    plt.ylabel('# reads')
    plt.xlabel('read length')
    plt.savefig(png_name)

    joined_pairs = hist1.read_count.sum()
    total_bp = (hist1.length*hist1.read_count).sum()
    mean_len = total_bp/joined_pairs
    
    if fastq2 is None:
        original_count = read_count(fastq1)/2
    else:
        original_count = read_count(fastq1)

    with open(outname, "w") as oh:
        print("metagenome", name, sep="\t", file=oh)
        print("seqs_before_join", original_count, sep="\t", file=oh)
        print("joined_pairs", joined_pairs, sep="\t", file=oh)
        print("total_bp", total_bp, sep="\t", file=oh)
        print("mean_length", mean_len, sep="\t", file=oh)
    return joined_pairs, original_count


def read_count(fname):
    """Count the number of reads and write metadata .count file.

    Args:
        fname (str): fastq or fasta file path

    Returns:
        read_count (int): integer of number of reads within fasta/fastq file
    """
    total = 0
    fq = True
    for name, seq, qual in readfx(fname):
        if not qual:
            fq = False
        break

    if fname.endswith("gz"):
        count_file = fname.rsplit(".gz", 1)[0] + ".count"
        cat = "gzip -d -c"
    else:
        count_file = fname + ".count"
        cat = "cat"

    if not fq:
        cmd = '%s %s | grep -c "^>"' % (cat, fname)
    else:
        cmd = '%s %s | wc -l' % (cat, fname)

    for line in run(cmd, description="Counting reads", iterable=True):
        total = int(line.rstrip())
        if fq:
            assert total % 4 == 0, "Multi-line or invalid FASTQ"
            total = int(total / 4)
    return total


def read_size_filter(fastx, readsize, outfile, cores=1):
    '''Read size filter

    Args:
        fastx (str): path to input fastq file
        readsize (int): minimum read size to keep
        outfile (str): location of outfile
    outputs:
        new fastq/fasta file with filtered reads

    '''
    if not outfile.endswith('.gz'):
        out = outfile
        outfile = outfile + '.gz'
    else:
        out=outfile.replace('.gz', "")

    if os.path.exists(outfile):
        return outfile

    with open(out, "w") as oh:
        totalreads = 0
        passedreads = 0
        for n, s, q in readfx(fastx):
            totalreads += 1
            if len(s) >= readsize:
                passedreads += 1
                if "fastq" in outfile:
                    print("@"+n, s, "+",q, sep="\n", file=oh)
                else:
                    print(">"+n, s, sep="\n", file=oh)
    print("{passedreads} out of {totalreads} passed the length filter.".format(**locals()))
    outfile = pigz_file(out, cores)
    return outfile



def compare_read_counts(joined_pairs, original_count):
    '''compare the number of joined pairs to the original number of reads in the forward fastq file

    Args:
        joined_pairs (int): the output of join_stats
        original_fastq (str): path to the original fastq file
    '''
    difference = original_count - joined_pairs
    if difference > original_count/2:
        print("ALERT! Joined library is less than half the size of original library.")
        #return False
    #else:
    #    return True
    return "there were {original_count} read pairs and {joined_pairs} joined reads".format(**locals())

def bbmap_stderr_todict(stderr):
    name = lambda i: i.split("\t")[0].strip().replace(":","").replace(" ","_")
    numreads = lambda i: i.split("\t")[2].strip()
    info = lambda i: i.split("\t")[1]
    ri = {}    # recruit info dict
    for i, l in enumerate(stderr.split("\n")):
        if len(l.split("\t")) == 5:
            ri[name(l)] = numreads(l)
        if len(l.split("\t")) == 2 or len(l.split("\t")) == 3:
            ri[name(l)] = info(l)
    return ri


def bbmap_stderr_info1(stderr, name):
    desired_info = ['perfect_best_site',
                    'unambiguous',
                    'Percent_scaffolds_with_any_coverage', 
                    'Percent_of_reference_bases_covered', 
                    'Average_coverage', 
                    'Reads_Used']
    indict = bbmap_stderr_todict(stderr)
    numreads = lambda i: i.split("\t")[2].strip()
    ser = Series([indict.get(i, "NA") for i in desired_info], index = desired_info, name=name)
    return ser

def bb_cov(ref, reads1, reads2=None, max_len=700, aln_id=0.95, threads=20, mem="20g", statsfile = True, outdir="", out_prefix=None):
    '''Run bbmap of designaged reads against reference fasta
    
    Args:
        ref (str): path to reference fasta
        reads1 (str): path to reads file
        reads2 (str): path to second reads file if data is paired, defaults to None
        aln_id (numeric between 0 and 1): percent identity to use to align reads, default is 0.95
        threads (int): how many threads to use for run, default is 20
        mem (str): how much memory to use in the syntax #unit, default is 20g which seems to work well for SAGs
        statsfile (boolean): if True, create coverage statistics file, if False, just return standard error; default is True
        outdir (str): designate if output data to be sent to directory other than current directory
        out_prefix (str): use to label output files.  Defaults to <reads>_bbmap2_<reference>

    Output:
        tuple of (bbmap standard error as text string, coverage statistics file)
    '''
    if statsfile:
        if out_prefix is None:
            rname = os.path.basename(ref).split(".")[0]
            mgname = os.path.basename(reads1).split(".")[0]
            out_prefix = "{mgname}_bbmap_{rname}".format(**locals())

        if reads2 is None:
            reads2cmd = ""
        else:
            reads2cmd = "in2={reads2}".format(**locals())
        
        oh = os.path.join(outdir, "{out_prefix}.covstats".format(**locals()))
        with file_transaction(oh) as tx:
            cmd = "bbmap.sh in={reads1} {reads2cmd} maxlen={max_len} minid={aln_id} covstats={tx} ref={ref} threads={threads} -Xmx{mem}".format(**locals())
            print(cmd)
            sterr = get_stderr(cmd)
        return sterr, oh
    else:
        cmd = "bbmap.sh in={reads1} {reads2cmd} minid={aln_id} ref={ref} threads={threads} -Xmx{mem}".format(**locals())
        print(cmd)
        sterr = get_stderr(cmd)
        return sterr

@cli.command('join', short_help='join meteganomic reads')
@click.argument('prefix')
@click.argument('outdir')
@click.argument('fq1', type=click.Path(exists=True))
@click.option('--fq2', default=None, help='input reverse reads file if reads are not interleaved')
@click.option('--mmd', type=click.FLOAT, default=0.05, help='mismatch density')
@click.option('--mino', type=click.INT, default=35, help='minimum overlap')
@click.option('--maxo', type=click.INT, default=150, help='maximum overlap')
@click.option('--threads', type=click.INT, default=20, help='numer of cores to run')
@click.option('--outdir', type=click.Path(), default="", help='output directory')
def join(prefix, fq1, fq2, mmd, mino, maxo, threads, outdir):
    '''Join metagenomic reads using flash
    Args:
        prefix (str): path with prefix included for location of output files
        fq1 (str): path to input fastq file, if interleaved, that's all you need, if separate forward and reverse include fq1
        options included in program description
    Output:
        joined read files, and statistics file of joined read process.
    '''
    outpath = os.path.join(outdir, prefix)
    outfiles = run_flash(outpath, fq1, fastq2=fq2, mismatch_density=mmd, min_overlap=mino,
        max_overlap=maxo, cores=threads)
    histin = outfiles[3]
    joined, reads = join_stats(histin, fq1, fastq2=fq2, prefix=prefix, outdir=outdir)
    compare_read_counts(joined, reads)


@cli.command('recruit', short_help='recruit metagenomic reads to template')
@click.option('--ref', type=click.Path(exists=True), default=None, help='sequence to recruit reads to in *.fasta format')
@click.option('--reads1', type=click.Path(exists=True), default=None, help='reads file in .fastq format')
@click.option('--reads2', type=click.Path(exists=True), default=None, help='second reads input if paired reads')
@click.option('--aln_id', type=click.FLOAT, default=0.95, help='minimum read alignment identity')
@click.option('--threads', type=click.INT, default=20, help='number of threads to use')
@click.option('--mem', default='20g', help='amount of memory to use')
@click.option('--outdir', type=click.Path(exists=True), default="./", help='output directory')
def recruit_with_stats(ref, reads1, reads2, aln_id, threads, mem, outdir):
    if read1 is None:
        raise IOError("Please provide reads file")
    if ref is None:
        raise IOError("Please provide reference sequence")
        
    sterr, cov_file = bbcov(ref, reads1, reads2=reads2, aln_id=aln_id, threads=threads, mem=mem, 
        statsfile=True, outdir=outdir, out_prefix=None)
    print(bbmap_stderr_info1(sterr))


@cli.command('recruit_multiref', short_help='recruit metagenomic reads to many references')
@click.argument('reads1', type=click.Path(exists=True))
@click.argument('refdir', type=click.Path(exists=True))
@click.option('--reads2', type=click.Path(exists=True), default=None, help='second reads input if paired reads')
@click.option('--aln_id', type=click.FLOAT, default=0.95, help='minimum read alignment identity')
@click.option('--threads', type=click.INT, default=20, help='number of threads to use')
@click.option('--mem', default='20g', help='amount of memory to use')
@click.option('--outdir', type=click.Path(), default="", help='output directory')
def recruit_multiref(reads1, refdir, reads2, aln_id, threads, mem, outdir):
    df = pd.DataFrame(index= ['perfect_best_site',
                'unambiguous',
                'Percent_scaffolds_with_any_coverage', 
                'Percent_of_reference_bases_covered', 
                'Average_coverage', 
                'Reads_Used'])
    
    if op.exists(outdir) == False:
        os.mkdir(outdir)

    readsname = os.path.basename(reads1).split(".")[0]
    out_table = os.path.join(outdir, readsname+"_recruit_info.txt")
    reflist = glob.glob(os.path.join(refdir+"*.fasta"))
    for ref in reflist:
        sterr, cov_file = bb_cov(ref, reads1, reads2=reads2, aln_id=aln_id, threads=threads, mem=mem, 
        statsfile=True, outdir=outdir, out_prefix=None)
        with open(os.path.join(outdir, "rawstderr.txt"),"a") as oh:
            print(sterr, file=oh)
        refname = re.sub(r"\.fa.*", "", os.path.basename(ref))
        ser = bbmap_stderr_info1(sterr, name=refname)
        df = pd.concat([df, ser], axis=1)
    df.to_csv(out_table, sep="\t")


@cli.command('jrmr', short_help='join reads then recruit to multiple references')
@click.option('--join', is_flag=True, help='if designated, join reads before recruitment')
@click.option('--len_filter', default=None, type=click.INT, help='minimum read length to keep for read recruitment')
@click.option('--max_len', default=600)
@click.option('--prefix', default=None, help='prefix to designate input metagenome, defaults to "mg"')
@click.option('--fq1', default=None, help='input reads in fastq form')
@click.option('--refdir', default=None, help='directory containing SAG genomes in *.fasta format')
@click.option('--fq2', default=None, help='input reverse reads file if reads are not interleaved')
@click.option('--mmd', type=click.FLOAT, default=0.05, help='mismatch density')
@click.option('--mino', type=click.INT, default=35, help='minimum overlap')
@click.option('--maxo', type=click.INT, default=150, help='maximum overlap')
@click.option('--threads', type=click.INT, default=20, help='number of threads to use')
@click.option('--aln_id', type=click.FLOAT, default=0.95, help='minimum read alignment identity')
@click.option('--mem', default='20g', help='amount of memory to use, must be in format <number><units>, e.g. "20g"')
@click.option('--outdir', default="", help='output directory')
def join_multirecruit(prefix, refdir, fq1, fq2, join, len_filter, max_len, mmd, mino, maxo, threads, outdir, aln_id, mem):
    if op.exists(outdir) == False:
        os.mkdir(outdir)
    
    print(fq1)

    if fq1 is None:
        raise IOError("Please designate input metagenome")

    if refdir is None:
        raise IOError("Please designate directory containing SAG reference fasta sequences")

    if len(glob.glob(os.path.join(refdir, "*.fasta"))) == 0:
        raise IOError("Reference directory contains no reference sequences.")

    if prefix is None:
        prefix = ".".join(os.path.basename(fq1).split(".")[:-1])
        print("prefix designator for input metagenome will be {prefix}".format(**locals()))
    
    #numbers = open(prefix+"_numbers.txt", "w")
    #perform join step if join designated
    if join:
        outpath = os.path.join(outdir, prefix)
        if os.path.exists(outpath+".extendedFrags.fastq.gz"):
            reads = outpath+".extendedFrags.fastq.gz"
            reads1 = reads
            reads2 = None
        else:
            outfiles = run_flash(outpath, fq1, fastq2=fq2, mismatch_density=mmd, min_overlap=mino,
                max_overlap=maxo, cores=threads)
            histin = outfiles[3]
            joined, reads = join_stats(histin, fq1, fastq2=fq2, prefix=prefix, outdir=outdir)
            print(compare_read_counts(joined, reads))
            
            reads1 = outfiles[0]
            reads2 = None
         
    else:
        reads1 = fq1
        reads2 = fq2

    # length filter reads if a number designated
    if len_filter is not None:
        if reads2 is not None:
            raise IOError("length filter not possible with forward and reverse reads")
        
        out_name = os.path.join(os.path.abspath(outdir), os.path.basename(reads1).replace(".fa", "_gt{len_filter}.fa".format(**locals())))
        
        if op.exists(out_name):
            reads1=out_name
        else:
            reads1 = read_size_filter(reads1, len_filter, out_name, cores=threads)

    #prep for bbcov:
    readsname = "_".join(os.path.basename(reads1).split(".")[:-1])
    out_table = os.path.join(outdir, readsname+"_recruit_info.txt")
    reflist = glob.glob(os.path.join(refdir+"*.fasta"))
    
    df = pd.DataFrame(index= ['perfect_best_site',
                'unambiguous',
                'Percent_scaffolds_with_any_coverage', 
                'Percent_of_reference_bases_covered', 
                'Average_coverage', 
                'Reads_Used'])   

    for ref in reflist:
        sterr, cov_file = bb_cov(ref=ref, reads1=reads1, reads2=reads2, aln_id=aln_id, threads=threads, mem=mem, 
        statsfile=True, outdir=outdir, out_prefix=None, max_len=max_len)
        with open(os.path.join(outdir, "rawstderr.txt"),"a") as oh:
            print(sterr, file=oh)
        refname = os.path.basename(ref).split(".")[0]
        ser = bbmap_stderr_info1(sterr, name=refname)
        df = pd.concat([df, ser], axis=1)
    
    df.to_csv(out_table, sep="\t")



if __name__ == '__main__':
    cli()



