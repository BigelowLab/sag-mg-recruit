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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scgc.utils import *
from scgc.fastx import readfx

'''
join reads using flash, filter reads less than length threshold
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
        list of output files in the order: joined reads, not combined fwd reads, not combined rev reads, hist, histogram
    """

    outsuffix = [".extendedFrags.fastq.gz", "out.notCombined_1.fastq.gz", ".notCombined_2.fastq.gz", ".hist", ".histogram"]
    outfiles = [prefix+i for i in outsuffix]
    for o in outfiles:
        if op.exists(o):
            exists = True
        else:
            exists = False
            break

    if exists == True:
        return outfiles

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
        passedreads = 0
        for n, s, q in readfx(outfile):
            passedreads += 1
        print("read filter output alreadyfound, {passedreads} are present in the filtered file".format(**locals()))
        return outfile, passedreads

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
    return outfile, passedreads



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


def join(prefix, fq1, fq2=None, mmd=.05, mino=35, maxo=350, threads=20, outdir=""):
    '''Join metagenomic reads using flash
    Args:
        prefix (str): path with prefix included for location of output files
        fq1 (str): path to input fastq file, if interleaved, that's all you need, if separate forward and reverse include fq2
        fq2 (str): path to referse read fastq file
        mmd (float<1): mismatch density for flash
        mino (int): minimum bp overlap
        maxo (int): maximum bp overlap
        threads (int): number of threads
        outdir (str): output directory path
    Output:
        joined read files, and statistics file of joined read process.
    '''
    outpath = os.path.join(outdir, prefix)
    outfiles = run_flash(outpath, fq1, fastq2=fq2, mismatch_density=mmd, min_overlap=mino,
        max_overlap=maxo, cores=threads)
    histin = outfiles[3]
    joined, reads = join_stats(histin, fq1, fastq2=fq2, prefix=prefix, outdir=outdir)
    compare_read_counts(joined, reads)
    return outfiles[0]


def process_multi_mgs(intable, outdir, threads, mmd, mino, maxo, minlen):
    '''Join and length filter a list of mg reads based on template from ../data/mg_template.csv
    
    Joins reads calling flash

    Args:
        intable (str): path to location of template file
        outdir (str): path to output directory
        mmd (float): mismatch density for flash (join program)
        mino (int): minimum bp overlap for joined reads
        maxo (int): maximum bp overlap
        threads (int): number of threads to run on
        outdir (str): desired directory for output files
        minlen (str): minimum read length to keep
    Returns:
        directory of result files, output table with number of reads and path to new input files added
    '''
    if op.exists(outdir) == False:
        safe_makedir(outdir)
    
    try:
        df = pd.read_csv(intable)
    except IOError as e:
        raise IOError("input table not found")

    df['name'] = [op.basename(i).split(".")[0] for i in df['mg_f']]
    mglist = list(df.mg_f[df['join']==False])
    
    # join reads identified as joined 
    tojoin = df.loc[df['join']==True]
    
    for n, f, r in zip(tojoin['name'], tojoin['mg_f'], tojoin['mg_r']):
        if pd.isnull(r):  # if reverse read cell is blank, but join=True, reads assumed to be interleaved
            r=None
        joinedfq = join(n, f, fq2=r, threads=threads, mmd=mmd, mino=mino, maxo=maxo, outdir=outdir)
        mglist.append(joinedfq)
    
    # size filter all reads:
    processed_mgs = []
    read_counts = []
    nameorder = []
    for m in mglist:
        nameorder.append(op.basename(m).split(".")[0])
        newfilename = op.basename(m).replace(".fastq", ".minlen_%s.fastq" % minlen)
        outfile = op.join(outdir, newfilename)
        new_out, count = read_size_filter(m, minlen, outfile, cores=threads)
        processed_mgs.append(new_out)
        read_counts.append(count)
    
    # create dataframe of results and merge with original info
    data = {'name':nameorder, 'to_recruit': processed_mgs, 'read_count':read_counts}
    newinfo = pd.DataFrame(data)
    res_table = pd.merge(df, newinfo, how='outer', on='name')
    tbl_name = op.join(outdir, "multi_mg_qc.csv")
    res_table.to_csv(tbl_name, sep=",")
    return res_table


@cli.command('join', short_help='join meteganomic reads')
@click.argument('prefix')
@click.argument('outdir')
@click.argument('fq1', type=click.Path(exists=True))
@click.option('--fq2', default=None, help='input reverse reads file if reads are not interleaved')
@click.option('--mmd', type=click.FLOAT, default=0.05, help='mismatch density')
@click.option('--mino', type=click.INT, default=35, help='minimum overlap')
@click.option('--maxo', type=click.INT, default=150, help='maximum overlap')
@click.option('--threads', type=click.INT, default=20, help='number of cores to run')
def run_join(prefix, outdir, fq1, fq2, mmd, mino, maxo, threads):
    join(prefix, fq1, fq2, mmd, mino, maxo, threads, outdir)


@cli.command('filt_size', short_help='remove reads below a certain size')
@click.option('--fastx', help='fasta or fastq file to filter')
@click.option('--outfile', default=None, help='name of output file with reads that passed the size threshold')
@click.option('--min_size', default=150, help='minimum read size to keep')
@click.option('--cores', default=1, help='number of cores to use')
def run_read_size_filt(fastx, outfile, min_size, cores):
    if outfile == None:
        outfile = ".".join(op.basename(fastx).split(".")[:-1])
    read_size_filter(fastx, readsize, outfile, cores=1)


@cli.command('process_table', short_help='process multiple metagenomes based on info supplied in input table')
@click.option('--table', type=click.Path(exists=True), help='input table, see wiki for format')
@click.option('--outdir', help='directory location to place output files')
@click.option('--threads', default=20, type=click.INT, help='number of cores to use')
@click.option('--mmd', type=click.FLOAT, default=0.05, help='mismatch density')
@click.option('--mino', type=click.INT, default=35, help='minimum overlap')
@click.option('--maxo', type=click.INT, default=150, help='maximum overlap')
@click.option('--minlen', type=click.INT, default=150, help='minimum read length')
def run_process_multi_mgs(table, outdir, threads, mmd, mino, maxo, minlen):
    process_multi_mgs(table, outdir, threads, mmd, mino, maxo, minlen)


if __name__ == '__main__':
    cli()

