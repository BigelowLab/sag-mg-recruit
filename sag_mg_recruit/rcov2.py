#!usr/bin/env python3
from __future__ import print_function

import contextlib
import subprocess
import os
import os.path as op
import pandas as pd
import itertools
import click
import logging
import pysam
import gzip
import contextlib
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
import functools
import tempfile
import shutil
from io import TextIOWrapper
from sarge import capture_stdout, capture_stderr



'''
append "real coverage" values to the end of the names of assembled contigs

requirements:
bwa
samtools
bedtools
'''

logger = logging.getLogger(__name__)

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """read recruitment."""
    pass

## from scgcpy dev branch for python 3:
def run(cmd, description=None, iterable=None, retcodes=[0]):
    """Runs a command using check_call.
    Args:
        cmd (str): shell command as a string
        description (Optional[str]):
        iterable (Optional[str]):
    """
    logging.debug("$> %s" % cmd)
    if description:
        logging.info(description)
    if iterable:
        return stdout_iter(cmd)
    else:
        p = capture_stderr(cmd)
        if p.returncode not in retcodes:
            for line in TextIOWrapper(p.stderr):
                logging.error(line.strip())
            raise subprocess.CalledProcessError(p.returncode, cmd=cmd)
        else:
            return p


@contextlib.contextmanager
def file_transaction(*rollback_files):
    """
    Wrap file generation in a transaction, moving to output if finishes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi", ".fastq.gz": ".count"}
    safe_names, orig_names = _flatten_plus_safe(rollback_files)
    # remove any half-finished transactions
    remove_files(safe_names)
    try:
        if len(safe_names) == 1:
            yield safe_names[0]
        else:
            yield tuple(safe_names)
    # failure -- delete any temporary files
    except:
        remove_files(safe_names)
        remove_tmpdirs(safe_names)
        raise
    # worked -- move the temporary files to permanent location
    else:
        for safe, orig in zip(safe_names, orig_names):
            if os.path.exists(safe):
                shutil.move(safe, orig)
                for check_ext, check_idx in exts.items():
                    if safe.endswith(check_ext):
                        safe_idx = safe + check_idx
                        if os.path.exists(safe_idx):
                            shutil.move(safe_idx, orig + check_idx)
        remove_tmpdirs(safe_names)


def remove_tmpdirs(fnames):
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)


def remove_files(fnames):
    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)


def _flatten_plus_safe(rollback_files):
    """
    Flatten names of files and create temporary file names.
    """
    tx_files, orig_files = [], []
    for fnames in rollback_files:
        if isinstance(fnames, str):
            fnames = [fnames]
        for fname in fnames:
            basedir = safe_makedir(os.path.dirname(fname))
            tmpdir = safe_makedir(tempfile.mkdtemp(dir=basedir))
            tx_file = os.path.join(tmpdir, os.path.basename(fname))
            tx_files.append(tx_file)
            orig_files.append(fname)
    return tx_files, orig_files


def safe_makedir(dname):
    """
    Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname


def file_exists(fnames):
    """
    Check if a file or files exist and are non-empty.
    parameters
        fnames : file path as string or paths as list; if list, all must exist
    returns
        boolean
    """
    if isinstance(fnames, str):
        fnames = [fnames]
    for f in fnames:
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            return False
    return True


@contextlib.contextmanager
def tmp_dir():
    d = None
    try:
        d = tempfile.mkdtemp()
        yield d
    finally:
        if d:
            shutil.rmtree(d)


##these functions are for bwa's illumina read aligner... not really applicable to the current project
def _aln(ref, fastq, tmp="/tmp", threads=8, threshold=0.05):
    sai = os.path.join(tmp, '%09d.sai' % random.randrange(0, 1e10))
    with file_transaction(sai) as tx:
        cmd = ("bwa aln -n {threshold} -t {threads} "
               "{ref} {fastq} > {tx}").format(**locals())
        run(cmd)
    return sai


def sampe_aln(ref, reads, bam_sorted, tmp="/tmp", threads=1, threshold=0.05):
    r1, r2 = tmp_split_reads(reads, tmp)
    with bwa_index(ref) as bwaidx:
        r1_sai = _aln(bwaidx, r1, tmp, threads, threshold)
        r2_sai = _aln(bwaidx, r2, tmp, threads, threshold)
        sampe = ("bwa sampe {ref} {r1_sai} {r2_sai} {r1} {r2} "
                 "| samtools view -bSF0x0004 - "
                 "| samtools sort -f -m 8 - {bam_sorted}").format(ref=bwaidx,
                                                       r1_sai=r1_sai,
                                                       r2_sai=r2_sai,
                                                       r1=r1,
                                                       r2=r2,
                                                       bam_sorted=bam_sorted)
        with file_transaction(bam_sorted) as tx:
            run(sampe)
        return bam_sorted

def samse_aln(ref, reads, bamout, tmp="/tmp", threads=8, threshold=0.05):
    
    with bwa_index(ref) as bwaidx:
        r_sai = _aln(bwaidx, reads, tmp, threads, threshold)
        samse = ("bwa samse {ref} {r_sai} {reads} | samtools view -bSF0x0004 - "
                 "| samtools sort -f -m 8 - {bam_sorted}")
        with file_transaction(bam_sorted) as tx:
            run(samse)
        return bam_sorted

## begin useful functions:
# from scgcpy.apps
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
        logger.info("creating bwa index for %s" % ref)
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
            logger.info("indexing bamfile %s" % bam_file)
    return bam_index


def _match_len(md):
    '''Calculate length of perfect alignment from md value

    Args:
        md (str): md string
    Returns:
        length (int): number of bp that were exact matches within alignment
    '''
    length=0
    number=""
    for i, c in enumerate(md):
        try:
            val = int(c)
            number = number + c
        except:
            if len(number) > 0:
                length += int(number)
                number=""
        if i == len(md)-1:
            length += int(number)
    return length


def filter_bam(bam, outbam, pctid = 95):
    '''Use sam "md" string to calculate pctid, filter aligned reads less than designated pctid
    
    Args:
        bam (str): path to bam file
        outbam (str): path to output bam file
        pctid (numeric between 0 and 100): percent identity desired for filtering
    Returns:
        tuple: path to (outbam, outfile)
    Outputs: 
        filtered bam file, output read count text file with format bamname \t # good reads
    '''
    with pysam.AlignmentFile(bam, "rb") as ih, pysam.AlignmentFile(outbam, "wb", template=ih) as oh:
        good = 0
        total = 0
        name = op.basename(outbam).split(".")[0]
        outfile = op.join(op.dirname(outbam),"%s.aln_count" % name)
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
            print(name, good, sep="\t", file=oh)
        logger.info("there were %s good read alignments out of %s total alignments" % (good, total))
    return outbam, outfile


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
        logger.info('running bwa mem: %s' % cmd)
        run(cmd)
        index_bam(tx_out_file)

    return out_file


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
        cmd = ("bedtools genomecov -dz -ibam {bam_file} > {tx_oh}").format(**locals())
        logger.info('running bedtools to get genome coverage: {cmd}'.format(cmd=cmd))
        subprocess.check_call(cmd, shell=True)
    return bedout
    

def cov_dict(bed_cov):
    '''
    create dict of mean coverage per contig based on bedtools coverage per base output
    '''
    coverage = pd.read_csv(bed_cov, sep="\t", header=None)
    return coverage.groupby([0])[2].mean().to_dict()
    
    
gzopen = lambda f: gzip.open(f) if f.endswith(".gz") else open(f)
 
 
def read_fasta(file_handle):
    for header, group in itertools.groupby(file_handle, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq
            

def append_rcov(infasta, outfasta, cov_dict):
    '''
    append "rcov" to the end of contig names within contig fasta file
    Args:
        infasta (path): input fasta file
        outfasta (path): output to write to
        cov_dict (dict): dictionary of new coverage values
    Returns:
        new fasta file with rcov appended to the end of the contig names
    '''
    with gzopen(infasta) as infile, open(outfasta, "w") as outfile:
        for name, seq in read_fasta(infile):
            cov = cov_dict.get(name, "NA")
            print(">%s_rcov_%.2f" % (name, cov), file = outfile)
            for i in xrange(0, len(seq), 80):        #wrap the lines
                print(seq[i:i + 80], file= outfile)


def bed_cov(fastq, reference, outdir, pctid, cores, cleanup, pe):
    fqpre = "_".join(op.basename(fastq).split(".")[:-1])
    ref_pre = "_".join(op.basename(reference).split(".")[:-1])
    outbam = op.join(os.path.abspath(outdir), fqpre+"_vs_"+ref_pre+".bam")
    
    if op.exists(outbam):
        bam = outbam
    elif pe:
        bam = bwa_mem(fastq, outbam, reference, options='-p', cores=cores)
    else:
        bam = bwa_mem(fastq, outbam, reference, options=None, cores=cores)             # run bwa mem 
    
    filtered_bam = bam.replace(".bam", "_{pctid}.bam".format(**locals()))
    countfile = filtered_bam.replace(".bam", ".alncount")

    if op.exists(filtered_bam) and op.exists(countfile):
        logger.warning('filtered bam file and countfile already exist, not re-running.')
        filtered_bam, goodcount = filtered_bam, countfile
    else:
        filtered_bam, goodcount = filter_bam(bam, bam.replace(".bam", "_{pctid}.bam".format(**locals())), pctid=pctid)

    bed = get_coverage(filtered_bam)                         # create per base coverage table
    logger.info("coverage table created, called: {bed}".format(bed=bed))
    
    if cleanup:
        #idx_files = [reference + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        #for f in idx_files+[bam, bam+".bai"]:
        for f in [bam, bam+".bai", filtered_bam]:
            logger.info('cleaning up extra bam files: {bam}, {bai}' 
                'and {fbam}'.format(bam=bam, bai=bam+'.bai', fbam=filtered_bam))
            os.remove(f)
    return bed, goodcount


def all_pairs(mglist, saglist):
    '''create list of tuples in which each mg is paired with every SAG'''
    allpairs = []
    for m in mglist:
        for s in saglist:
            allpairs.append((m, s))
    return allpairs


def split_threads_by_runs(threads, runs):
    '''given the number of threads designated to the program and the number of 
    runs to be done simultaneously, 
    return a pool process for threads designated 
    and the number of cores to use per run
    '''
    cores = threads//runs
    if cores == 0:
        raise IOError("number of threads must be greater than or equal to the number of runs")
    return cores


def bed_cov_it(fqref, outdir, pctid, cores, cleanup, pe=False):
    '''This function takes a fastq file and a reference, entered in the first input as a tuple, 
    1. aligns them using bwa
    2. filters them based on pctid alignment threshold
    3. creates a counts file of how many reads passed the filter
    4. creates a genome coverage table using bedtools
    5. returns a tuple of the coverage file and the count file
    Args:
        fqref: tuple as (path to metagenome fastq, path to reference SAG)
        outdir (str): output directory
        pctid (int): percent identity of aligned reads to keep
        cores (int): # threads to use 
        cleanup (boolean): if true, delete output bam and bai file after running
        pe (boolean): if true, reads are paired 
    Output:
        bedtools coverage file, alignment count file
    '''
    fastq = fqref[0]
    reference = fqref[1]
    fqpre = "_".join(op.basename(fastq).split(".")[:-1])
    ref_pre = "_".join(op.basename(reference).split(".")[:-1])
    outbam = op.join(os.path.abspath(outdir), fqpre+"_vs_"+ref_pre+".bam")
    
    if pe:
        bam = bwa_mem(fastq, outbam, reference, options='-p', cores=cores)
    else:
        bam = bwa_mem(fastq, outbam, reference, options=None, cores=cores)             # run bwa mem 
    bam, goodcount = filter_bam(bam, bam.replace(".bam", "_{pctid}.bam".format(**locals())), pctid=pctid)   # filter aligned reads based on pctid
    
    bed = get_coverage(bam)                         # create per base coverage table
    print("coverage_table_created, called:", bed)
    if cleanup:
        #idx_files = [reference + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        #for f in idx_files+[bam, bam+".bai"]:
        for f in [bam, bam+".bai"]:
            os.remove(f)
    return bed


def run_bed_cov_mp(fqref_list, outdir, pctid, total_cores, cleanup, pe):
    '''create bedtools coverage and counts file for each mg, sag pair designated in input list

    Args:
        fqref_list (list of tuples): list of tuples of (path to metagenome fastq, path to reference SAG)
        outdir (str): path to directory to place output files
        pctid (int): percent identity of aligned reads to keep
        cores (int): # threads to use 
        cleanup (boolean): if true, delete output bam and bai file after running
        pe (boolean): if true, reads are paired 
    Returns:
        list of output bedtools coverage files
    Outputs:
        two files in designated directory:
        .genomecoverage file (bedtools genome coverage) 
        .alncount file that reports the number of reads recruited to reference genome 
    '''
    
    core_pr = split_threads_by_runs(total_cores, 3)       # set up pool and number of cores per run
    #rbc = functools.partial(bed_cov_it, outdir=outdir, pctid=pctid, 
    #              cores=core_pr, cleanup=True)            # set up partial function to be fed into pool.map
    #p = multiprocessing.Pool(processes=total_cores)
    
    #triplets=[]
    #for i in range(0, len(fqref_list), 3):
    #    trip=[]
    #    for j in fqref_list[i:i+3]:
    #        trip.append(j)
    #    triplets.append(trip)

    #results = []
    #for tri in triplets:                   # run three alignments at a time 
    #    results += p.map(rbc, tri)
    #p.close()
    results=[]
    for fr in fqref_list:
        res = bed_cov(fr[0], fr[1], outdir=outdir, pctid=pctid, cores=total_cores, cleanup=True, pe=False)
        results.append(res)
    return results


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
    cols = ['sag',
            'metagenome',
            'Percent_scaffolds_with_any_coverage', 
            'Percent_of_reference_bases_covered', 
            'Average_coverage', 
            'total_reads_recruited']
    data = [sag, metagenome, 
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
    
@cli.command("append-coverage", short_help='append real coverage values to contig names')
@click.option('--fastq', help='input fastq file')
@click.option('--reference', help='reference contigs')
@click.option('--outfasta', help='name of output fasta file')
@click.option('--pe', is_flag=True, help='enter if fastq consists of interleaved reads')
@click.option('--cores', default=8, help='number of cores for bwa to use')
@click.option('--cleanup', default=False, help='boolean to designate if extra files created during run should be deleted')
def append_real_cov(fastq, reference, outfasta, cleanup, cores, pe):
    '''
    append real coverage values based on bwa alignment
    Args:
        fastq (path): reads file
        contigs (path): contigs fasta file
        outfasta (path): name of new contigs fasta file
        cleanup (boolean): if True, remove all sam/bam/bed files after creating new fasta
    Returns:
        new fasta file with _rcov_<mean_coverage> appended to the end of the seq names
    '''
    fqpre = "_".join(op.basename(fastq).split(".")[:-1])
    ref_pre = "_".join(op.basename(fastq).split(".")[:-1])
    outbam = op.join(os.path.abspath(outdir), fqpre+"_vs_"+ref_pre+".bam")
    
    if pe:
        bam = bwa_mem(fastq, outbam, reference, options='-p', cores=cores)
    else:
        bam = bwa_mem(fastq, outbam, reference, cores=cores) 

    bed = get_coverage(bam)      
    cd = cov_dict(bed)                               # dict of mean coverage per contig 
    append_rcov(reference, outfasta, cd)               # append real mean coverage (rcov) to contig name
    if cleanup:                                     
        idx_files = [reference + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        for f in idx_files+[bam, bam+".bai"]:
            if op.exists(f):
                os.remove(f)
    print("done.")


@cli.command("print-coverage", short_help='output bedtools coverage table')
@click.option('--fastq', help="input fastq file")
@click.option('--reference', help='reference contigs')
@click.option('--outdir', default="", help='output directory')
@click.option('--pe', is_flag=True, help='enter if fastq consists of interleaved reads')
@click.option('--pctid', default=95, help="minimum percent identity to keep")
@click.option('--cores', default=8, help='number of cores for bwa to use')
@click.option('--cleanup', default=False, help='boolean to designate if extra files created during run should be deleted') 
def run_print_real_cov(fastq, reference, outdir, pctid, cores, cleanup, pe):
    bed_cov(fastq, reference, outdir, pctid, cores, cleanup, pe)


@cli.command('multi-recruit', short_help='run bwa on all designated metagenomes against all designated SAGs')
@click.option('--fqlist_file', help='path to file containing a list of metagenome files in fastq format')
@click.option('--saglist_file', help='path to file containing a list of sag file in fasta format')
@click.option('--outdir', help='output directory')
@click.option('--pctid', default=95, help='minimum percent identity to keep for read alignment')
@click.option('--total_cores', default=12, help='total cores to use to run this script')
@click.option('--log', default="multirecruit.log", help='name of logfile for multi-recruit')
def run_multi_recruit(fqlist_file, saglist_file, outdir, pctid, total_cores, log):
    if op.exists(outdir)==False:
        safe_makedir(outdir)

    if log == None:
        log = op.join(outdir, "mulitrecruit.log")
        print("logfile is: %s" % log)

    logging.basicConfig(filename=log, level=logging.DEBUG)
    mglist = [i for i in open(fqlist_file).read().split("\n") if len(i) > 0]
    saglist = [i for i in open(saglist_file).read().split("\n") if len(i) > 0]
    outfile = op.join(outdir, "combined_recruitment_info.tsv")
    pairs = all_pairs(mglist, saglist)
    covtbls = run_bed_cov_mp(pairs, outdir, pctid, total_cores, cleanup=True, pe=False)
    combined = genome_cov_table(covtbls)
    combined.to_csv(outfile, sep="\t")
    print("result table written to {outfile}".format(**locals()))


if __name__ == '__main__':
    cli()

