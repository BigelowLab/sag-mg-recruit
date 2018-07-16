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
from itertools import groupby
import click
import os.path as op
import logging
import shutil
import gzip
import pysam
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distutils.spawn import find_executable

from scgc.utils import file_transaction, safe_makedir, run, tmp_dir, pigz_file


__version_info__ = (0, 0, 1)
__version__ = '.'.join(map(str, __version_info__))
REQUIRES = ["bedtools", "samtools", "checkm", "bwa", "gzip", "gunzip"]


logger = logging.getLogger(__name__)
# wgs_factors = {'illumina':0.8376, 'pyro':1}
gzopen = lambda i: gzip.open(i) if i.endswith(".gz") else open(i)


def check_dependencies(executables):
    exes = []
    for exe in executables:
        if not find_executable(exe):
            exes.append(exe)
    if len(exes) > 0:
        for exe in exes:
            print("`%s` not found in PATH." % exe)
        sys.exit(1)


def readfx(fastx):
    if not file_exists(fastx):
        logger.critical("File Not Found: %s" % fastx)
        raise IOError(2, "No such file:", fastx)

    fx = ""
    try:
        fx = pysam.FastxFile(fastx)
        for f in fx:
            yield f.name, f.sequence, f.quality
    finally:
        if fx:
            fx.close()


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

    outsuffix = [".extendedFrags.fastq.gz", ".notCombined_1.fastq.gz", ".notCombined_2.fastq.gz", ".hist", ".histogram"]
    outfiles = [prefix+i for i in outsuffix]
    for o in [outfiles[0], outfiles[3], outfiles[4]]:
        if op.exists(o):
            exists = True
        else:
            exists = False
            logger.debug("file %s does not exist" % o)
            break

    if exists:
        return outfiles

    print("FASTQ1 for FLASH is:", fastq1)
    with file_transaction(outfiles) as tx_outfiles:
        if fastq2 is None:
            logger.info("only one fastq file provided, assuming it is interleaved")
            cmd = ("flash --interleaved-input -x {mismatch} -m {min} "
                   "-M {max} -t {threads} -z -o {prefix} {fastq1}").format(mismatch=mismatch_density,
                        min=min_overlap, max=max_overlap, threads=cores, prefix=prefix,
                        fastq1=fastq1)
        else:
            cmd = ("flash -m {min} -M {max} -x {mismatch} -t {threads} "
           "-z -o {prefix} {fastq1} {fastq2}".format(mismatch=mismatch_density,
                        min=min_overlap, max=max_overlap, threads=cores,
                        fastq1=fastq1, fastq2=fastq2, prefix=prefix))
        logger.info("running flash: %s" % cmd)
        try:
            run(cmd, description="Joining reads with flash")
        except:
            logger.error("join step could not be performed for {fastq1}".format(fastq1=fastq1))
            return ["","","","","",""]
    # delete uncombined reads to save space
    for f in outfiles:
        if "notCombined" in f:
            os.remove(f)
    return outfiles


def join_stats(inhist, fastq1, fastq2=None, prefix="", outdir=""):
    '''
    Print join stats and png of read size distribution

    Args:
        inhist (str): path to .hist output created by flash
        fastq1
        fastq2
        prefix (str): location to write output to, defaults to current directory
        outdir

    Returns:
        list of [int, int] which are # joined pairs, original paired read count
    '''
    name = os.path.basename(inhist).replace(".hist", "")

    outname = os.path.join(outdir, "{prefix}_joinstats.txt".format(prefix=prefix))
    png_name = os.path.join(outdir, "{prefix}_joinstats.png".format(prefix=prefix))

    hist1 = pd.read_table(inhist, header=None)
    hist1.columns = ["length", "read_count"]
    fig = plt.plot(hist1['length'], hist1['read_count'], color='b')
    plt.ylabel('# reads')
    plt.xlabel('read length')
    plt.savefig(png_name)

    joined_pairs = hist1.read_count.sum()
    total_bp = (hist1.length * hist1.read_count).sum()
    mean_len = total_bp / joined_pairs

    if fastq2 is None:
        original_count = read_count(fastq1, outdir)[0] / 2
    else:
        original_count = read_count(fastq1, outdir)[0]

    with open(outname, "w") as oh:
        print("metagenome", name, sep="\t", file=oh)
        print("seqs_before_join", original_count, sep="\t", file=oh)
        print("joined_pairs", joined_pairs, sep="\t", file=oh)
        print("total_bp", total_bp, sep="\t", file=oh)
        print("mean_length", mean_len, sep="\t", file=oh)
    return joined_pairs, original_count


def read_count(fname, directory, minlen=0):
    """Count the number of reads and write metadata .count file.

    Args:
        fname (str): fastq or fasta file path

    Returns:
        read_count (int): integer of number of reads within fasta/fastq file
    """
    total_reads = 0
    total_bp = 0
    fq = True
    for name, seq, qual in readfx(fname):
        if not qual:
            fq = False
        break
    if op.exists(fname) == False:
        logger.error("could not find file: %s" % fname)
        return 0

    if fname.endswith("gz"):
        count_file = op.join(directory, "{0}_minlen{1}.count".format(".".join(op.basename(fname).split(".")[:-2]), minlen))
    else:
        count_file = op.join(directory, "{0}_minlen{1}.count".format(".".join(op.basename(fname).split(".")[:-1]), minlen))

    if op.exists(count_file):
        total_reads, total_bp = open(count_file).read().split("\n")[0:2]
        return total_reads, total_bp

    for name, seq, qual in readfx(fname):
        if len(seq) >= minlen:
            total_reads += 1
            total_bp += len(seq)
        else:
            continue

    with open(count_file, "w") as oh:
        print(total_reads, total_bp, file=oh, sep="\n")
    return total_reads, total_bp


def read_size_filter(fastx, read_size, outfile, cores=1):
    '''Read size filter

    Args:
        fastx (str): path to input fastq file
        read_size (int): minimum read size to keep
        outfile (str): location of outfile

    Returns:
        List of [str, int]: FASTQ file path of output and total passing reads.
    '''
    if not outfile.endswith('.gz'):
        out = outfile
        outfile = outfile + '.gz'
    else:
        out=outfile.replace('.gz', "")

    if os.path.exists(outfile):
        passed_reads = 0
        for n, s, q in readfx(outfile):
            passed_reads += 1
        # log
        logger.info("read filter output already found, {passed_reads} are present in the filtered file".format(passed_reads=passed_reads))
        return outfile, passed_reads

    with open(out, "w") as oh:
        total_reads = 0
        passed_reads = 0
        for n, s, q in readfx(fastx):
            total_reads += 1
            if len(s) >= read_size:
                passed_reads += 1
                if "fastq" in outfile:
                    print("@"+n, s, "+",q, sep="\n", file=oh)
                else:
                    print(">"+n, s, sep="\n", file=oh)
    logger.info(("for {fastx}, {passed_reads} out of {total_reads} passed the length filter. "
                 "printing to {outfile}").format(fastx=fastx,
                                                 passed_reads=passed_reads,
                                                 total_reads=total_reads,
                                                 outfile=outfile))
    outfile = pigz_file(out, cores)
    return outfile, passed_reads


def compare_read_counts(joined_pairs, original_count):
    '''compare the number of joined pairs to the original number of reads in the forward fastq file

    Args:
        joined_pairs (int): the output of join_stats
        original_fastq (str): path to the original fastq file

    Returns:
        str
    '''
    difference = int(original_count) - joined_pairs
    if difference > int(original_count) / 2:
        logger.warning("ALERT! Joined library is less than half the size of original library.")
    else:
        logger.debug("number of joined reads is greater than half the size of the original library")
    return "there were {original_count} read pairs and {joined_pairs} joined reads".format(original_count=original_count, joined_pairs=joined_pairs)


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

    Returns:
        path to joined reads (str)
    '''
    outpath = os.path.join(outdir, prefix)
    outfiles = run_flash(outpath, fq1, fastq2=fq2, mismatch_density=mmd, min_overlap=mino,
        max_overlap=maxo, cores=threads)
    if len(outfiles[0]) == 0:
        return outfiles[0]
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
        pandas.DataFrame - output table with number of reads
    '''
    if op.exists(op.join(outdir, "multi_mg_qc_minlen{}.txt".format(minlen))):
        logger.info("looks like these metagenomes have already been processed, loading mg stats table.")
        logger.info("If the metagenomes have been moved, please delete {} and rerun sag-mg-recruit".format(op.join(outdir, "multi_mg_qc_minlen{}.txt".format(minlen))))
        return pd.read_csv(op.join(outdir, "multi_mg_qc_minlen{}.txt".format(minlen)))

    if op.exists(outdir) == False:
        safe_makedir(outdir)


    df = pd.read_csv(intable)

    to_recruit = []
    for i, r in df.iterrows():
        if r.join:
            # there's a mix of op and os.path calls
            to_recruit.append(op.join(outdir, "%s.extendedFrags.fastq.gz" % r['name']))
        else:
            to_recruit.append(r.mg_f)
    df['to_recruit'] = to_recruit

    # join reads identified as joined
    # there's a mix of to_recruit and tojoin without the underscore. def go with the underscore!
    to_join = df.loc[df['join']==True]

    for n, f, r in zip(to_join['name'], to_join['mg_f'], to_join['mg_r']):
        # if reverse read cell is blank, but join=True, reads assumed to be interleaved
        if pd.isnull(r):
            r = None

        joined_fq = join(n, f, fq2=r, threads=threads, mmd=mmd, mino=mino, maxo=maxo, outdir=outdir)
        logger.info("reads joined for {}".format(n))

    total_counts = [read_count(m, outdir, 0) for m in df['to_recruit']]
    df['total_reads_before_filter'] = [c[0] for c in total_counts]
    df['total_bp_before_filter'] = [c[1] for c in total_counts]

    len_counts = [read_count(m, outdir, minlen) for m in df['to_recruit']]
    df['read_count'] = [c[0] for c in len_counts]
    df['bp_count'] = [c[1] for c in len_counts]
    # or
    # df['read_count'] = df['to_recruit'].apply(read_count) or maybe it's df['to_recruit'].apply(lambda x: read_count(x))
    # create dataframe of results

    df.to_csv(op.join(outdir, "multi_mg_qc_minlen{}.txt".format(minlen)), sep="\t", index=False)
    return df


# SAG functions
def mask_sag(input_gb, out_fasta):
    '''output masked contigs with rRNA changed to 'N's

    Args:
        input_gb (str): path to annotated input genbank formatted genome
        out_fasta (str): where to write the output fasta to

    Returns:
        str - fasta file with rRNA regions masked with 'N's
    '''
    #if input_gb.endswith(".gb") == False or input_gb.endswith(".gbk") == False:
    #    logger.error("input file does not appear to be in genbank format.  Please check.")
    #    return None
    if op.exists(out_fasta):
        return out_fasta
        logger.info("masked fasta for {} already exists".format(input_gb))
    with open(input_gb, "rU") as input_handle, open(out_fasta, "w") as oh:
        rrna_count = 0
        for r in SeqIO.parse(input_handle, "genbank"):
            print(">", r.name, sep="", file=oh)
            s = r.seq
            cloc = []
            masked = ""
            for f in r.features:
                if f.type == "rRNA" or f.type == "RNA":
                    if ('gene' in f.qualifiers and
                            ("16S" in str(f.qualifiers['gene']).upper() or
                            "23S" in str(f.qualifiers['gene']).upper() or
                            "Subunit Ribosomal RNA".upper() in str(f.qualifiers['gene']).upper() or
                            "suRNA".upper() in str(f.qualifiers['gene']).upper())):
                        # if the 'type' is rRNA, it should be masked... don't have to check for 16 or 23S
                        cloc.append(f.location)
                        logger.info('rRNA gene found on contig %s' % r.name)
                        rrna_count += 1
                    elif ('product' in f.qualifiers and
                            ("16S" in str(f.qualifiers['product']).upper() or
                            "23S" in str(f.qualifiers['product']).upper() or
                            "Subunit Ribosomal RNA".upper() in str(f.qualifiers['product']).upper() or
                            "suRNA".upper() in str(f.qualifiers['product']).upper())):
                        # if the 'type' is rRNA, it should be masked... don't have to check for 16 or 23S
                        cloc.append(f.location)
                        logger.info('rRNA gene found on contig %s' % r.name)
                        rrna_count += 1
                elif f.type == "misc_feature":
                    if ('note' in f.qualifiers and
                        ("16S" in str(f.qualifiers['note']).upper() or
                         "23S" in str(f.qualifiers['note']).upper())):
                        cloc.append(f.location)
                        logger.info('miscellaneous feature with "16S" or "23S" in note found')
                        rrna_count += 1

            # if the contig contains one rRNA gene (most common if rRNA present)
            if len(cloc) == 1:
                logger.debug("contig {name} has 1 rRNA gene".format(name=r.name))
                masked += s[0:cloc[0].start - 1]
                masked += 'N'*(cloc[0].end - cloc[0].start)
                masked += s[cloc[0].end - 1:]

            elif len(cloc) > 1:
                logger.debug("contig {name} has more {num} rRNA genes".format(name=r.name, num=len(cloc)))
                for i in range(0, len(cloc)):
                    # if it's the first entry
                    if i == 0:
                        masked += s[0:cloc[i].start - 1]
                        masked += 'N'*(cloc[i].end - cloc[i].start)
                    # if it's the last entry
                    elif i == len(cloc) - 1:
                        masked += s[cloc[i - 1].end - 1:cloc[i].start - 1]
                        masked += 'N'*(cloc[i].end - cloc[i].start)
                        masked += s[cloc[i].end:]
                    else:
                        masked += s[cloc[i - 1].end -1 :cloc[i].start - 1]
                        masked += 'N'*(cloc[i].end - cloc[i].start)
            # if no rRNA on contig, just return the sequence, unmasked
            else:
                logger.debug("contig {name} does not have any annotated rRNA genes".format(name=r.name))
                masked = s

            for i in range(0, len(masked), 80):
                print(masked[i:i + 80], file=oh)
    logger.info('%s rRNA genes found in %s' % (rrna_count, op.basename(input_gb)))
    return out_fasta


def gbk_to_fasta(input_gb, out_fasta):
    with open(input_gb, "rU") as input_handle, open(out_fasta, "w") as oh:
        for r in SeqIO.parse(input_handle, "genbank"):
            print(">", r.name, sep="", file=oh)
            for i in range(0, len(r.seq), 60):
                print(r.seq[i:i+60], file=oh)
    return out_fasta


def process_gb_sags(tbl, outdir):
    '''process SAGs according to intsructions in input table

    Args:
        tbl (str): path to input table with the following columns:
        outdir (str): path to output directory.  does not have to exists before running function.

    Returns:
        list of re-named SAGs in fasta format, with 16/23S masked by 'N's if mask == True
    '''
    safe_makedir(outdir)

    fas_sags = []
    df = pd.read_csv(tbl)
    for i, l in df.iterrows():
        if l['mask'] == True:
            if l.gbk_file is not None and op.exists(l.gbk_file):
                outfasta = op.join(outdir, l.sag_name + ".masked.fasta")
                fas_sags.append(mask_sag(l.gbk_file, outfasta))
                logger.info("{} masked".format(l.sag_name))
            else:
                logger.error("could not find input genbank file to mask")
        # if mask not designated, write sag to fasta if gbk supplied, else use supplied fasta
        elif l['mask'] == False:
            out_fasta = op.join(outdir, l.sag_name + ".fasta")
            if l.fasta_file is None:
                fas_sags.append(gbk_to_fasta(l.gbk_file, out_fasta))
            else:
                with open(l.fasta_file) as ih, open(out_fasta, "w") as oh:
                    for name, seq in read_fasta(ih):
                        print(">{}".format(name), file=oh)
                        seq_fix = seq.replace(" ","")
                        for j in range(0, len(seq_fix), 60):
                            print(seq_fix[j:j+60], file=oh)
            fas_sags.append(out_fasta)
    return fas_sags


def sag_checkm_completeness(fasta, cores):
    '''run checkm lineage_wf on SAG to get completeness values

    Args:
        fasta (str): full path to SAG genomic contigs in fasta format
        cores (int): number of cores to use to run checkm

    Returns:
        "completeness" statistics as a pandas dataframe
    '''
    logger.info("Running checkm on %s " % fasta)

    fasta = op.abspath(fasta)

    with tmp_dir() as tdir:
        bindir = op.join(tdir, "bindir")
        safe_makedir(bindir)
        outdir = op.join(tdir, "outdir")
        safe_makedir(outdir)

        tmp_fasta = op.join(bindir, op.basename(fasta))

        shutil.copy(fasta, tmp_fasta)
        assert op.exists(tmp_fasta)
        logger.debug("{} created".format(tmp_fasta))

        completeness_tsv = op.join(outdir, "completeness.tsv")

        cmd = "checkm lineage_wf -f {outfile} --tab_table -q -x fasta -t {cores} {binpath} {outdir}".format(outfile=completeness_tsv, outdir=outdir, cores=cores, binpath=bindir)

        logger.info("running checkm lineage, command is: {cmd}".format(cmd=cmd))
        run(cmd)
        completeness = pd.read_csv(completeness_tsv, sep="\t", header=0)
    return completeness


def process_sag_fastas(saglist, outfile, cores, checkm):
    '''calculate checkM completeness value given the SAGs listed in a file

    Args:
        sagfile (str): path to file containing a list of paths, one per line, of SAG fasta files to analyze
        outfile (str): path to location where output table will be written
        cores (int): number of cores to use to run checkm
        checkm (bool): option for program.  If user opts out of running checkm, just return table with checkm
    Returns:
        tab-delimited file of checkm completeness per SAG
    '''
    logger.info("gathering checkM completeness values for all SAGs listed in file: {}".format(saglist))
    df = pd.DataFrame(columns=['Bin Id', 'Marker lineage', '# genomes', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity', 'total_bp'])

    if checkm == "True" or checkm == True:
        logger.info("checkm completeness calculations requested")
        df = pd.DataFrame(columns=['Bin Id', 'Marker lineage', '# genomes', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity', 'total_bp'])
        for s in saglist:
            if not op.exists:
                logger.error("SAG not found for %s" % s)
                continue
            if not op.isfile:
                logger.error("%s is not a file" % s)
                continue
            completeness = sag_checkm_completeness(s, cores=cores)
            if completeness is None:
                logger.info("completeness stats for %s not determined" % s)
                continue

            length = count_fasta_bp(s)
            logger.debug("sag %s is %s bp in length" % (s, length))

            completeness['total_bp'] = length
            # completeness['calculated_length'] = int(completeness.total_bp * 100/completeness.Completeness)

            df = pd.concat([df, completeness])
    else:
        binid = []
        completeness = []
        total_bp = []
        for s in saglist:
            if not op.exists:
                logger.error("SAG not found for %s" % s)
                continue
            if not op.isfile:
                logger.error("%s is not a file" % s)
                continue
            binid.append(op.basename(s))
            completeness.append("NA")
            total_bp.append(count_fasta_bp(s))
        df = pd.DataFrame(data={'Bin Id':binid, 'Completeness':completeness, 'total_bp':total_bp})

    df.to_csv(outfile, sep="\t", index=False)
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


def count_fasta_bp(sagfasta):
    total_length = 0
    with gzopen(sagfasta) as infile:
        for name, seq in read_fasta(infile):
            total_length += len(seq)
    return total_length


## coverage functions ##
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
    """Build an index for a bam file.

    Args:
        bam_file (str): alignment file path

    Returns:
        string: index file path
    """
    bam_index = bam_file + '.bai'
    if not file_exists(bam_index):
        with file_transaction(bam_index) as tx_out_file:
            run('samtools index %s %s' % (bam_file, tx_out_file))
    return bam_index


def read_overlap_pctid(l, pctid, min_len, overlap=0):
    real_len = l.infer_query_length()
    aln_len = l.query_alignment_length
    mismatch = l.get_tag("NM")

    aln_overlap = (aln_len / real_len) * 100
    aln_pctid = ((aln_len - mismatch) / aln_len) * 100
    if aln_overlap >= overlap and aln_pctid >= pctid and aln_len >= min_len:
        return True
    else:
        return False


def filter_bam(bam, outbam, pctid=95, minlen=150, overlap=0,):
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as ih, pysam.AlignmentFile(outbam, "wb", template=ih) as oh:
        good = 0
        good_bp = 0
        total = 0
        name = op.basename(outbam).split(".")[0]
        outfile = ".".join(outbam.split(".")[:-1]) + ".aln_count"
        for i, l in enumerate(ih):
            if l.is_duplicate:
                continue

            total += 1
            #md = l.get_tag("MD")
            #match = _match_len(md)
            #pct_match = (match)/l.rlen * 100

            #if pct_match > pctid:
            #    good += 1
            #    oh.write(l)
            if read_overlap_pctid(l, pctid, minlen, overlap):
                good += 1
                good_bp += l.query_alignment_length
                oh.write(l)

        with open(outfile, "w") as oh:
            print(good, good_bp, sep="\n", file=oh)
        logger.info("for %s, there were %s good read alignments out of %s total alignments" % (bam, good, total))
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
    # basestring is python 2 specific
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
        reference : path to bwa index
        options : bwa mem options
        cores : int
    returns
        output file path : string
    """

    if file_exists(out_file):
        return out_file

    assert op.exists(fastq), "Could not find fastq file {}".format(fastq)
    assert op.exists(reference), "Could not find reference file {}".format(reference)

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
               "-ShuF4q2 -@ {cores} - | samtools sort -m 8G -@ {cores} -o {result}"
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
    # basestring is python 2 specific
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
        cmd = ("bedtools genomecov -d -ibam {bam} > {out}").format(bam=bam_file, out=tx_oh)
        subprocess.check_call(cmd, shell=True)
    return bedout


def print_real_cov(fastq, reference, outdir, pctid, overlap, minlen, cores, cleanup, mgname=None, referencename=None, pe=None):
    ''' calculate per-base coverage using bwa, samtools and bedtools

    Args:
        fastq (path): input metagenome in fastq format
        reference (path): input reference genome in fasta format
        outdir (path): path to output directory
        pctid (int): minimum percent identity for aligned reads
        cores (int): number of cores to run on
        cleanup (boolean): if True, delete .bam and .bai files after coverage is calculated
        pe: intput True if reads are paird and interleaved
    Outputs:
        if cleanup = False, outputs .bam, .bai, .genomecoverage, and an .aln_count file
        if cleanup = True, outputs .genomecoverage and .aln_count file only
    Returns:
        path to genome coverage file
    '''
    if mgname is None:
        fqpre = op.basename(fastq).split(".")[0]
    else:
        fqpre = mgname

    if referencename is None:
        ref_pre = op.basename(reference).split(".")[0]
    else:
        ref_pre = referencename

    outbam = op.join(os.path.abspath(outdir), fqpre + "_vs_" + ref_pre + ".bam")
    # print(outbam)

    if pe:
        bam = bwa_mem(fastq, outbam, reference, options='-p', cores=cores)
    else:
        bam = bwa_mem(fastq, outbam, reference, options=None, cores=cores)

    bam = filter_bam(bam,
                     bam.replace(".bam", ".pctid{pctid}.overlap{overlap}.minlen{minlen}.bam".format(pctid=pctid,
                        overlap=overlap, minlen=minlen)), overlap=overlap, pctid=pctid, minlen=minlen)
    # create per base coverage table
    bed = get_coverage(bam)
    logger.info("coverage_table_created, called:{}".format(bed))

    if cleanup:
        idx_files = [reference + x for x in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        for f in idx_files + [bam, bam + ".bai"]:
            if op.exists(f):
                os.remove(f)
    return bed


def get_recruit_info(gcov):
    '''calculate information on recruited reads based on bedtools genomecoverage table

    Args:
        gcov (str): path to genome coverage file with recruitment pipeline naming convention of:
            metagenome_vs_sag.genomecoverage
            metagenome_vs_sag.aln_count file must also exists within the same directory

    Returns:
        pandas dataframe of genome coverage statistics
    '''
    countfile = gcov.replace("genomecoverage", "aln_count")
    with open(countfile) as infile:
        recruit_count, recruit_bp = infile.read().split("\n")[0:2]
        #recruit_count = infile.read().split()[1].strip()

    metagenome = op.basename(gcov).split("_vs_")[0]
    sag = op.basename(gcov).split("_vs_")[1].split(".")[0]
    cols = ['sag',
            'metagenome',
            'Percent_scaffolds_with_any_coverage',
            'Percent_of_reference_bases_covered',
            'Average_coverage',
            'total_reads_recruited',
            'total_bp_recruited']

    try:
        coverage = pd.read_csv(gcov, sep="\t", header=None)
    except:
        logger.warning("no genome coverage data for {sag}-{metagenome} recruitment".format(sag=sag, metagenome=metagenome))
        data = [sag, metagenome, 0, 0, 0, 0, 0]
        return pd.DataFrame(data, index=cols).transpose()

    mean_per_contig = coverage.groupby([0])[2].mean()
    sum_per_contig = coverage.groupby([0])[2].sum()
    contig_size = coverage.groupby([0])[1].max() + 1
    mean_sag_coverage = mean_per_contig.mean()
    totalbp = contig_size.sum()

    uncovered_bp = len(coverage[coverage[2] == 0])
    pct_covered = ((totalbp - uncovered_bp) / totalbp) * 100
    total_scaffold = len(sum_per_contig)
    uncovered_contig = len(sum_per_contig[sum_per_contig == 0])
    pct_scaffolds_covered = ((total_scaffold - uncovered_contig) / total_scaffold) * 100

    data = [sag,
           metagenome,
           pct_scaffolds_covered,
           pct_covered,
           mean_sag_coverage,
           recruit_count,
           recruit_bp]
    return pd.DataFrame(data, index=cols).transpose()



def genome_cov_table(gcov_listm dev=False):
    '''create large dataframe of metagenome recruitment information given a number of genome coverage files

    Args:
        gcov_list (list): list of paths to genome coverage files

    Returns:
        pandas dataframe summary of recruitment information from all genome coverage files in list
    '''

    big = pd.concat([get_recruit_info(g) for g in gcov_list])
    if not dev:
        for g in gcov_list:
            os.remove(g)
    return big


def cov_from_list(fastqlist, referencelist, mg_names, reference_names, outdir, pctid, overlap, minlen, cores, outtable, cleanup=False, keep_coverage=False):
    '''create large datframe of recruitment information for multiple SAGs against multiple metagenomes
    Args:
        fastqlist(list): list of paths to metagenomic reads in fastq format
        referencelist(list): list of paths to SAGs in fasta format
        outdir (str): name of output directory to print genome coverage to
        pctid (int): percent identity of read alignments to keep
        cores (int): number of cores to use
        outtable (str): path to result table
        cleanup (boolean): if true, delete all bam-type files
    Returns:
        writes output table to outtable
        a pandas dataframe of the result table
    '''
    bedlist = []

    if mg_names is None:
        mg_names = [op.basename(i).split(".")[0] for i in fastqlist]

    if reference_names is None:
        reference_names = [op.basename(i).split(".")[0] for i in referencelist]

    for fn, f in zip(mg_names, fastqlist):
        for rn, r in zip(reference_names, referencelist):
            bed = print_real_cov(f, r, outdir=outdir, pctid=pctid, overlap=overlap, minlen=minlen, cores=cores, cleanup=cleanup, mgname=fn, referencename=rn)
            bedlist.append(bed)

    table = genome_cov_table(bedlist, dev=keep_coverage)
    table.to_csv(outtable, sep="\t", index=False)
    logger.info("result table written to {outfile}".format(outfile=outtable))
    return table

def concatenate_fastas(fastalist, outfasta):
    with open(outfasta, "w") as oh:
        for s in fastalist:
            with gzopen(s) as ih:
                for name, seq in read_fasta(ih):
                    print(">" + name, file=oh)
                    # wrapping the sequence in fasta output
                    for i in range(0, len(seq), 80):
                        print(seq[i:i + 80], file=oh)
    return outfasta


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('input_mg_table')
@click.argument('input_sag_table')
# global options
@click.option('--outdir',
              default=None,
              help='directory location to place output files')
@click.option('--cores',
              default=8,
              show_default=True,
              help='number of cores to run on')
# mg processing options
@click.option('--mmd',
              type=click.FLOAT,
              default=0.05,
              show_default=True,
              help='for join step: mismatch density')
@click.option('--mino',
              type=click.INT,
              default=35,
              show_default=True,
              help='for join step: minimum overlap')
@click.option('--maxo',
              type=click.INT,
              default=150,
              show_default=True,
              help='for join step: maximum overlap')
# coverage options
@click.option('--minlen',
              type=click.INT,
              default=150,
              show_default=True,
              help='for alignment and mg read count: minimum alignment length to include; minimum read size to include')
@click.option('--pctid',
              default=95,
              show_default=True,
              help="for alignment: minimum percent identity to keep within overlapping region")
@click.option('--overlap',
              default=0,
              show_default=True,
              help="for alignment: percent read that must overlap with reference sequence to keep")
@click.option('--log',
              default=None,
              help='name of log file, else, log sent to standard out')
@click.option('--concatenate',
               type=click.BOOL,
               default=True,
               show_default=True,
               help='include concatenated SAG in analysis')
@click.option('--checkm',
                type=click.BOOL,
                default=True,
                show_default=True,
                help='should checkm be run on the SAGs?')
@click.option('--keep_coverage',
                type=click.flag,
                help='if you want to keep the genome coverage table (large)')
def main(input_mg_table, input_sag_table, outdir, cores,
         mmd, mino, maxo, minlen, pctid, overlap, log, concatenate, checkm, keep_coverage=False):
    if log is None:
        log = logging.StreamHandler(sys.stdout)
        log.setLevel(logging.INFO)
        logger.addHandler(log)
    else:
        logging.basicConfig(filename=log, level=logging.INFO)

    check_dependencies(REQUIRES)
    parms = str(print("PARAMETERS for sag-mg-recruit:",
                  "input_mg_table = {}".format(input_mg_table),
                  "input_sag_table = {}".format(input_sag_table),
                  "outdir = {}".format(outdir),
                  "cores = {}".format(cores),
                  "mmd = {}".format(mmd),
                  "mino = {}".format(mino),
                  "maxo = {}".format(maxo),
                  "minlen = {}".format(minlen),
                  "pctid = {}".format(pctid),
                  "overlap = {}".format(overlap),
                  "log = {}".format(log),
                  "concatenate = {}".format(concatenate),
                  "checkm={}".format(checkm), sep="\n"))

    logger.info(parms)

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
    summaryout = op.join(outdir, "summary_table_pctid{pctid}_minlen{minlen}_overlap{overlap}.txt".format(pctid=pctid, minlen=minlen, overlap=overlap))

    logger.info("processing the metagenomes")
    tbl_name = op.join(mgdir, "multi_mg_qc_minlen{}.txt".format(minlen))
    if op.exists(tbl_name):
        mgtbl = pd.read_csv(tbl_name, sep="\t")
        logger.info("Metagenomes have already been processed.  Loading {}".format(tbl_name))
    else:
        mgtbl = process_multi_mgs(input_mg_table, mgdir, threads=cores, mmd=mmd, mino=mino, maxo=maxo, minlen=minlen)

    mglist = mgtbl['to_recruit']
    mgnames = mgtbl['name']

    logger.info("processing sag table")
    saglist = process_gb_sags(input_sag_table, sagdir)
    logger.info("calculating SAG completeness using CheckM")

    completeness_out = op.join(sagdir, "sag_completeness.txt")

    if op.exists(completeness_out):
        sagtbl = pd.read_csv(completeness_out, sep="\t")
        logger.info("SAGs have already been processed.  Loading {}".format(completeness_out))
    else:
        sagtbl = process_sag_fastas(saglist, completeness_out, cores, checkm)

    logger.info("running bwa read recruitment")

    if concatenate == "True" or concatenate == True:
        sagconcat = op.join(sagdir, "concatenated_sags.fasta")
        if op.exists(sagconcat) == False:
            sagconcat = concatenate_fastas(saglist, sagconcat)

        saglist.append(sagconcat)

    coverage_out = op.join(covdir, "coverage_info_pctid{pctid}_minlen{minlen}_overlap{overlap}.txt".format(pctid=pctid, minlen=minlen, overlap=overlap))
    if op.exists(coverage_out):
        covtbl = pd.read_csv(coverage_out, sep="\t")
        logger.info("bwa recruitment has already been done. loading {}".format(coverage_out))
    else:
        covtbl = cov_from_list(mglist, saglist, mgnames, None, covdir, pctid, overlap, minlen, cores, coverage_out, cleanup=False, keep_coverage=keep_coverage)

    # process tables to make summary table
    logger.info('putting together summary tables')
    sagtbl['sag'] = [i.split(".")[0] for i in sagtbl['Bin Id']]
    sagtbl.rename(columns={'Completeness': 'sag_completeness',
                           'total_bp': 'sag_total_bp'},
                  inplace=True)
    sagshort = sagtbl[['sag', 'sag_completeness', 'sag_total_bp']]

    mgtbl.rename(columns={'name': 'metagenome',
                          'wgs_technology': 'mg_wgs_technology',
                          'read_count': 'mg_read_count', 'bp_count':'mg_bp_count'},
                 inplace=True)
    mgshort = mgtbl[['metagenome', 'mg_wgs_technology', 'mg_read_count', 'mg_bp_count']]

    covtbl['sag'] = [i.split(".")[0] for i in covtbl['sag']]

    summary = covtbl.merge(mgshort, how='outer', on='metagenome')
    summary = summary.merge(sagshort, how='outer', on='sag')

    summary[['sag_total_bp', 'total_reads_recruited', 'total_bp_recruited', 'mg_read_count', 'mg_bp_count']] = summary[['sag_total_bp', 'total_reads_recruited', 'total_bp_recruited', 'mg_read_count', 'mg_bp_count']].convert_objects(convert_numeric=True)


    summary['sag_size_mbp'] = summary.sag_total_bp / 1000000
    summary['MG_reads_per_SAG_mbp'] = summary.total_reads_recruited / summary.sag_size_mbp
    summary['MG_bp_per_SAG_mpb'] = summary.total_bp_recruited / summary.sag_size_mbp
    summary['prop_total_MG_reads_recruited_per_SAG_mbp'] = summary.MG_reads_per_SAG_mbp / summary.mg_read_count
    summary['prop_total_MG_bp_recruited_per_SAG_mbp'] = summary.MG_bp_per_SAG_mpb / summary.mg_bp_count


    summary.to_csv(summaryout, sep="\t", index=False)

    logger.info('process completed.')
    return summary


if __name__=='__main__':
    main()
