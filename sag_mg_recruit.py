'''
Determine representation of SAGs within metagenomes 

read recruitment of metagenomes to SAGs using bwa.
'''

from __future__ import print_function
from __future__ import division

import subprocess
import os
import sys
import pandas as pd
from pandas import Series
from itertools import groupby
import click
import glob
import re
import os.path as op
import logging
import shutil
import gzip
import pysam
from Bio import SeqIO
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distutils.spawn import find_executable


from scgc.utils import *
from scgc.fastx import readfx


__version_info__ = (0, 0, 0)
__version__ = '.'.join(map(str, __version_info__))
REQUIRES = ["bedtools", "samtools", "checkm", "bwa", "gzip", "gunzip", "Rscript"]


logger = logging.getLogger(__name__)
wgs_factors = {'illumina':0.8376, 'pyro':1}


def check_dependencies(executables):
    exes = []
    for exe in executables:
        if not find_executable(exe):
            exes.append(exe)
    if len(exes) > 0:
        for exe in exes:
            print("`%s` not found in PATH." % exe)
        sys.exit(1)


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
        else:
            cmd = ("flash -m {min} -M {max} -x {mismatch} -t {threads} "
           "-z -o {prefix} {fastq1} {fastq2}".format(mismatch=mismatch_density,
                        min=min_overlap, max=max_overlap, threads=cores,
                        fastq1=fastq1, fastq2=fastq2, prefix=prefix))
        logger.info("running flash: %s" % cmd)
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
        print("read filter output already found, {passedreads} are present in the filtered file".format(**locals()))
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
    logger.info("for {fastx}, {passedreads} out of {totalreads} passed the length filter."
                "printing to {outfile}".format(**locals()))
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
        logger.warning("ALERT! Joined library is less than half the size of original library.")
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
        wd = os.getcwd()
        f = f.replace(wd, "./")
        r = r.replace(wd, "./")
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
    res_table.to_csv(tbl_name, sep="\t")
    return res_table


# SAG functions
def sag_checkm_completeness(fasta,  cores):
    '''run checkm lineage_wf on SAG to get completeness values

    Args:
        fasta (str): full path to SAG genomic contigs in fasta format
        cores (int): number of cores to use to run checkm
    Returns:
        "completeness" statistics as a pandas dataframe
    '''
    logger.info("Running checkm on %s " % fasta)

    fasta = op.abspath(fasta)
    if op.isdir == True or op.exists == False:
        return None
    
    with tmp_dir() as tdir:
        bindir = op.join(tdir, "bindir")
        safe_makedir(bindir)
        outdir = op.join(tdir, "outdir")
        safe_makedir(outdir)
        
        tmp_fasta = op.join(bindir, op.basename(fasta))
        
        try:
            shutil.copy(fasta, tmp_fasta)
            assert op.exists(tmp_fasta)
            print(tmp_fasta, "created")
        except Exception, e:
            print("copying %s to the temporary directory failed, %s" % (fasta, e))
            return None

        logger.info("Running lineage workflow on %s" % fasta)
        
        completeness_tsv = op.join(outdir, "completeness.tsv")
        
        cmd = "checkm lineage_wf -f {outfile} --tab_table -q -x fasta -t {cores} {binpath} {outdir}".format(outfile=completeness_tsv, outdir=outdir, cores=cores, binpath=bindir)
        
        logger.info("running checkm lineage, command is: {cmd}".format(**locals()))
        run(cmd)
        completeness = pd.read_csv(completeness_tsv, sep="\t", header=0)
    return completeness
    


def checkm_completeness(saglist, outfile, cores):
    '''calculate checkM completeness value given the SAGs listed in a file

    Args:
        sagfile (str): path to file containing a list of paths, one per line, of SAG fasta files to analyze
        outfile (str): path to location where output table will be written
    Returns:
        tab-delimited file of checkm completeness per SAG
    '''
    logger.info("gathering checkM completeness values for all SAGs listed in file: saglist")
    df = pd.DataFrame(columns=['Bin Id', 'Marker lineage', '# genomes', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity', 'total_bp'])

    for s in saglist:
        if op.exists == False:
            logger.error("SAG not found for %s" % s)
            continue
        if op.isfile == False:
            logger.error("%s is not a file" % s)
            continue
        completeness = sag_checkm_completeness(s, cores=cores)
        if completeness is None:
            logger.info("completeness stats for %s not determined" % s)
            #print("completeness stats for %s not determined" % s)
            continue

        length = count_fasta_bp(s)
        #print("sag %s is %s bp in length" % (s, length))

        completeness['total_bp'] = length
        completeness['calculated_length'] = int(completeness.total_bp * 100/completeness.Completeness)
        
        df = pd.concat([df, completeness])

    df.to_csv(outfile, sep="\t")
    return df


def read_fasta(file_handle):
    '''Fasta iterator'''
    for header, group in groupby(file_handle, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


gzopen = lambda i: gzip.open(i) if i.endswith(".gz") else open(i)


def count_fasta_bp(sagfasta):
    total_length = 0
    with gzopen(sagfasta) as infile:
        for name, seq in read_fasta(infile):
            total_length += len(seq)
    return total_length



## coverage functions
def bwa_index(reference):
    """Builds an index using `bwa index`.
    Args:
        reference (str): file path of reference fasta
    Returns:
        str: file path of reference as it's used as the prefix in `bwa index`
    """
    ref = op.abspath(reference)
    idx_files = [ref + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
    if not file_exists(idx_files):
        cmd = "bwa index %s" % ref
        subprocess.check_call(cmd, shell=True)
    return reference


def index_bam(bam_file):
    """
    Build an index for a bam file.
    parameters
        bam_file : alignment file path
    returns
        index file name : string
    """
    bam_index = bam_file + '.bai'
    if not file_exists(bam_index):
        with file_transaction(bam_index) as tx_out_file:
            run('samtools index %s %s' % (bam_file, tx_out_file))
    return bam_index


def _match_len(md):
    length=0
    mismatch=0
    number=""
    for i, c in enumerate(md):
        try:
            val = int(c)
            number = number+c
        except:
            if len(number) > 0:
                length += int(number)
                number=""
        if i == len(md)-1:
            length += int(number)
    return length


def filter_bam(bam, outbam, pctid = 95):
    with pysam.AlignmentFile(bam, "rb") as ih, pysam.AlignmentFile(outbam, "wb", template=ih) as oh:
        good = 0
        total = 0
        name = op.basename(outbam).split(".")[0]
        outfile = ".".join(outbam.split(".")[:-1])+".aln_count"
        for i, l in enumerate(ih):
            if l.is_duplicate:
                continue

            total += 1
            md = str(l).split("\t")[-1].split(",")[3].replace(")","").replace("'","").strip()  # get md value from raw bam entry
            match = _match_len(md)
            pct_match = (match)/l.rlen * 100

            if pct_match > pctid:
                good += 1
                oh.write(l)
        with open(outfile, "w") as oh:
            print(name, good, file=oh)
        print("there were %s good read alignments out of %s total alignments" % (good, total))
    return outbam


def _remove_option(options, item, flag=False):
    """
    remove item from options. item is assumed to not be a flag, therefore
    two sequential args are removed from options.
    parameters
        options : list of args
        item : string to search for among options
        flag : whether for not item is boolean
    returns
        options : list
    """
    if item in options:
        x = options.index(item)
        options.pop(x)
        if not flag:
            options.pop(x)
    return options


def filter_options(options, predefined_options):
    """
    filter predefined options from user-specified options.
    parameters
        options : string or list of options for application
        predefined_options : list of tuples specifying option and whether or
            not the option is a boolean flag
    returns
        options : list
    """
    options = options.split() if isinstance(options, basestring) else options
    for option, flag in predefined_options:
        options = _remove_option(options, option, flag)
    return options


def bwa_mem(fastq, out_file, reference, options, cores=1):
    """
    align reads using bwa mem.
    parameters
        fastq : path to reads
        out_file : path to aligned reads bam
        index : path to bwa index
        options : bwa mem options
        cores : int
    returns
        output file path : string
    """
    if file_exists(out_file):
        return out_file
    predefined_options = [('-t', False)]
    
    if options is not None:
        options = filter_options(options, predefined_options)
        opts = " ".join(options)
    else:
        opts = ""
    
    logger.info("Mapping %s to %s using bwa mem" % (fastq, reference))
    
    reference = bwa_index(reference)

    with file_transaction(out_file) as tx_out_file:
        cmd = ("bwa mem -t {cores} {options} {index} {fastq} | samtools view "
               "-ShuF4q2 - | samtools sort -o -m 8G - tmp > {result}"
              ).format(cores=cores,
                       options=opts,
                       index=reference,
                       fastq=fastq,
                       result=tx_out_file)
        run(cmd)
        index_bam(tx_out_file)

    return out_file


def file_exists(fnames):
    """
    Check if a file or files exist and are non-empty.
    parameters
        fnames : file path as string or paths as list; if list, all must exist
    returns
        boolean
    """
    if isinstance(fnames, basestring):
        fnames = [fnames]
    for f in fnames:
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            return False
    return True


def get_coverage(bam_file, bedout=None):
    '''
    create per base coverage patterns from sorted bam
    ''' 
    bedgraph = ""
    filename, ext = op.splitext(bam_file)
    if bedout is None:
        bedout = filename + ".genomecoverage"
    
    if op.exists(bedout):
        return bedout

    with file_transaction(bedout) as tx_oh:
        cmd = ("bedtools genomecov -d -ibam {bam_file} > {tx_oh}").format(**locals())
        subprocess.check_call(cmd, shell=True)
    return bedout

### calculate coverage:

def print_real_cov(fastq, reference, outdir, pctid, cores, cleanup, pe):
    fqpre = op.basename(fastq).split(".")[0]
    ref_pre = op.basename(reference).split(".")[0]
    outbam = op.join(os.path.abspath(outdir), fqpre+"_vs_"+ref_pre+".bam")
    
    if pe:
        bam = bwa_mem(fastq, outbam, reference, options='-p', cores=cores)
    else:
        bam = bwa_mem(fastq, outbam, reference, options=None, cores=cores)             # run bwa mem 

    bam = filter_bam(bam, bam.replace(".bam", ".{pctid}.bam".format(**locals())), pctid=pctid)

    bed = get_coverage(bam)                         # create per base coverage table
    print("coverage_table_created, called:", bed)
    
    if cleanup:
        idx_files = [reference + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        for f in idx_files+[bam, bam+".bai"]:
            if op.exists(f):
                os.remove(f)
    return bed


def get_recruit_info(gcov):
    countfile = gcov.replace("genomecoverage", "aln_count")
    with open(countfile) as infile:
        recruit_count = infile.read().split()[1].strip()
        
    metagenome = op.basename(gcov).split("_vs_")[0]
    sag = op.basename(gcov).split("_vs_")[1].split(".")[0]
    try:
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
    except:
        mean_per_contig = 0
        sum_per_contig = 0
        contig_size = 0
        mean_sag_coverage = 0
        totalbp = 0
        uncovered_bp = 0
        pct_covered = 0
        total_scaffold = 0
        uncovered_contig = "NA"
        pct_scaffolds_covered = "NA" if uncovered_contig == 'NA' else (total_scaffold - uncovered_contig)/total_scaffold *100

    cols = ['sag',
            'metagenome',
            'Percent_scaffolds_with_any_coverage', 
            'Percent_of_reference_bases_covered', 
            'Average_coverage', 
            'total_reads_recruited']
    data = [sag, 
           metagenome, 
           pct_scaffolds_covered,
           pct_covered, 
           mean_sag_coverage,
           recruit_count]
    df = pd.DataFrame(data, index=cols).transpose()
    return df


def genome_cov_table(gcov_list):
    cols = ['sag',
            'metagenome',
            'Percent_scaffolds_with_any_coverage', 
            'Percent_of_reference_bases_covered', 
            'Average_coverage', 
            'total_reads_recruited']
    big = pd.DataFrame(columns=cols)
    for g in gcov_list:
        new=get_recruit_info(g)
        big = pd.concat([big, new])
    return big


def cov_from_list(fastqlist, referencelist, outdir, pctid, cores, outtable, cleanup=False, pe=None):
    bedlist = []
    for f in fastqlist:
        for r in referencelist:
            bed = print_real_cov(f, r, outdir=outdir, pctid=pctid, cores=cores, cleanup=cleanup, pe=pe)
            bedlist.append(bed)
    table = genome_cov_table(bedlist)
    table.to_csv(outtable, sep="\t")
    return table
    print("result table written to {outfile}".format(outfile=outtable))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('input_mg_table')
@click.argument('input_sag_list')
# global options
@click.option('--outdir', default=None, help='directory location to place output files')
@click.option('--cores', default=8, help='number of cores to run on')
# mg processing options
@click.option('--mmd', 
              type=click.FLOAT, 
              default=0.05, 
              help='for join step: mismatch density for join step in mg processing')
@click.option('--mino', 
              type=click.INT, 
              default=35, 
              help='for join step: minimum overlap for join step')
@click.option('--maxo', 
              type=click.INT, 
              default=150, 
              help='for join step: maximum overlap for join step')
@click.option('--minlen', 
              type=click.INT, 
              default=150, 
              help='for metagenomes: minimum read length')
# coverage options
@click.option('--pctid', 
              default=95, 
              help="for alignment: minimum percent identity to keep")
@click.option('--log', 
              default=None, 
              help='name of log file for this process')
def main(input_mg_table, input_sag_list, outdir, cores, 
         mmd, mino, maxo, minlen, pctid, log):
    if log == None:
        log = "sag-mg-recruit.log"
        print("Logfile is: %s" % log)
    logging.basicConfig(filename=log, level=logging.INFO)
    
    check_dependencies(REQUIRES)

    if outdir is None:
        outdir = "mg_sag_recruitment"
        logger.info("output directory is {outdir}".format(**locals()))
    outdir = safe_makedir(outdir)
    mgdir = op.join(outdir, 'mgs')
    mgdir = safe_makedir(mgdir)
    sagdir = op.join(outdir, 'sags')
    sagdir = safe_makedir(sagdir)
    covdir = op.join(outdir, 'coverage')
    covdir = safe_makedir(covdir)
    summaryout = op.join(outdir, "summary_table.txt")
    
    logger.info("processing the metagenomes")
    tbl_name = op.join(mgdir, "multi_mg_qc.csv")
    if op.exists(tbl_name):
        mgtbl = pd.read_csv(tbl_name, sep="\t")
        logger.info("Metagenomes have already been processed.  Loading {tbl_name}".format(**locals()))
    else:
        mgtbl = process_multi_mgs(input_mg_table, mgdir, threads=cores, mmd=mmd, mino=mino, maxo=maxo, minlen=minlen)
    
    mglist = mgtbl['to_recruit']
    
    logger.info("calculating SAG completeness using CheckM")
    saglist = [i for i in open(input_sag_list).read().split("\n") if len(i) > 0]
    completeness_out = op.join(sagdir, "sag_completeness.txt")
    
    if op.exists(completeness_out):
        sagtbl = pd.read_csv(completeness_out, sep="\t")
        logger.info("SAGs have already been processed.  Loading {completeness_out}".format(**locals()))
    else:
        sagtbl = checkm_completeness(saglist, completeness_out, cores)
    
    logger.info("running bwa read recruitment")
    coverage_out = op.join(covdir, "coverage_info.txt")
    if op.exists(coverage_out):
        covtbl = pd.read_csv(coverage_out, sep="\t")
        logger.info("bwa recruitment has already been done. loading {coverage_out}".format(**locals()))
    else:
        covtbl = cov_from_list(mglist, saglist, covdir, pctid, cores, coverage_out, cleanup=True)

    # process tables to make summary table
    logger.info('putting together summary tables')
    #sagtbl['sag'] = [i.split("_")[0] for i in sagtbl['Bin Id']]
    sagtbl['sag'] = [i.split(".")[0] for i in sagtbl['Bin Id']]
    sagtbl.rename(columns={'Completeness':'sag_completeness', 
                           'total_bp':'sag_total_bp', 
                           'calculated_length':'sag_calculated_length'}, 
                           inplace=True)
    sagshort = sagtbl[['sag', 'sag_completeness', 'sag_total_bp', 'sag_calculated_length']]

    mgtbl.rename(columns={'name':'metagenome', 
                          'wgs_technology':'mg_wgs_technology', 
                          'read_count':'mg_read_count'}, 
                          inplace=True)
    mgshort = mgtbl[['metagenome', 'mg_wgs_technology', 'mg_read_count']]
    
    covtbl['sag'] = [i.split(".")[0] for i in covtbl['sag']]

    summary = covtbl.merge(mgshort, how='outer', on='metagenome')
    summary = summary.merge(sagshort, how='outer', on='sag')    
    
    try:
        summary['prop_mg_recruited'] = (summary.total_reads_recruited*100/summary.sag_completeness)/summary.mg_read_count
        summary['prop_mg_adjusted'] = summary['prop_mg_recruited']*summary['mg_wgs_technology'].map(wgs_factors)
    except Exception as inst:
        logger.warning(type(inst))     # the exception instance
        logger.warnning(inst.args)      # arguments stored in .args
        logger.warning(inst)     
        logger.warning("prop_mg_recruited and prop_mg_adjusted unable to be calculated.")
        summary['prop_mg_recruited'] = "NA"
        summary['prop_mg_adjusted'] = "NA"
    
    summary.to_csv(summaryout, sep="\t")

    logger.info('process completed.')
    return summary


if __name__=='__main__':
    main()

