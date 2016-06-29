#!usr/bin/env python
from __future__ import print_function
from __future__ import division

import subprocess
import os
import os.path as op
import pandas as pd
import itertools
import click
import logging
import pysam
from Bio import SeqIO
import gzip
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scgc.utils import run, file_transaction, safe_makedir


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
        cmd = ("bedtools genomecov -dz -ibam {bam_file} > {tx_oh}").format(**locals())
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


@cli.command("append_cov", short_help='append real coverage values to contig names')
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
            os.remove(f)
    print("done.")
    


def print_real_cov(fastq, reference, outdir, pctid, cores, cleanup, pe):
    fqpre = "_".join(op.basename(fastq).split(".")[:-1])
    ref_pre = "_".join(op.basename(reference).split(".")[:-1])
    outbam = op.join(os.path.abspath(outdir), fqpre+"_vs_"+ref_pre+".bam")
    
    if pe:
        bam = bwa_mem(fastq, outbam, reference, options='-p', cores=cores)
    else:
        bam = bwa_mem(fastq, outbam, reference, options=None, cores=cores)             # run bwa mem 

    bam = filter_bam(bam, bam.replace(".bam", "_{pctid}.bam".format(**locals())), pctid=pctid)

    bed = get_coverage(bam)                         # create per base coverage table
    print("coverage_table_created, called:", bed)
    
    if cleanup:
        idx_files = [reference + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        for f in idx_files+[bam, bam+".bai"]:
            os.remove(f)
    return bed


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
    sag = op.basename(gcov).split("_vs_")[1].replace(".genomecoverage","").split("_")[0]
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


def cov_from_list(fastqlist, referencelist, outdir, pctid, cores, outtable, cleanup=False, pe=None):
    bedlist = []
    for f in fastqlist:
        for r in referencelist:
            bed = print_real_cov(f, r, outdir=outdir, pctid=pctid, cores=cores, cleanup=cleanup, pe=pe)
            bedlist.append(bed)
    table = genome_cov_table(bedlist)
    table.to_csv(outtable, sep="\t")
    return outtable
    print("result table written to {outfile}".format(outfile=outtable))


@cli.command("print_cov", short_help='output bedtools coverage table')
@click.option('--fastq', help="input fastq file")
@click.option('--reference', help='reference contigs')
@click.option('--outdir', default="", help='output directory')
@click.option('--pe', is_flag=True, help='enter if fastq consists of interleaved reads')
@click.option('--pctid', default=95, help="minimum percent identity to keep")
@click.option('--cores', default=8, help='number of cores for bwa to use')
@click.option('--cleanup', default=False, help='boolean to designate if extra files created during run should be deleted') 
def run_print_real_cov(fastq, reference, outdir, pctid, cores, cleanup, pe):
    print_real_cov(fastq, reference, outdir, pctid, cores, cleanup, pe)


@cli.command("print_cov_list", short_help='output recruitment stats for all' 
              'combinations of fastq and reference fastas provided in list files')
@click.option('--fastqlist', help='file containing list of fastq files, one per line')
@click.option('--reflist', help='file containing list of reference fastas, one per line')
@click.option('--outdir', default="", help='output directory')
@click.option('--pe', is_flag=True, help='enter if fastq consists of interleaved reads')
@click.option('--pctid', default=95, help="minimum percent identity to keep")
@click.option('--cores', default=8, help='number of cores for bwa to use')
@click.option('--cleanup', default=False, help='boolean to designate if extra files created during run should be deleted') 
@click.option('--outtable', default=None, help='name of output table to be used')
def run_cov_from_list(fastqlist, reflist, outdir, pe, pctid, outtable, cores, cleanup):
    if outtable is None:
        outtable = op.join(outdir, "recruitment_info.txt")
    else:
        outtable = op.join(outdir, outtable)
    mglist = [i for i in open(fastqlist).read().split("\n") if len(i) > 0]
    saglist = [i for i in open(reflist).read().split("\n") if len(i) > 0]
    outtable = cov_from_list(mglist, saglist, outdir, pctid, cores, outtable, cleanup, pe)
    print('coverage info written to {out}'.format(out=outtable))


if __name__ == '__main__':
    cli()