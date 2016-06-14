from __future__ import print_function

import os.path as op
from scgc.utils import tmp_dir, safe_makedir, run, file_exists
import logging
import pandas as pd
import shutil
import click
from itertools import groupby
import gzip
# import subprocess
# from sarge import get_stderr


logger = logging.getLogger(__name__)

'''
sag qc information for read recruitment
'''

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """completeness calculation via checkm."""
    pass


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
    if op.isfile == False or op.exists == False:
        return None
    #completeness = 0.0
    #number_unique_markers = 0
    #number_multi_copy = 0
    #taxonomy_contained = "NA"
    #taxonomy_sister_lineage = "NA"

    # we want the fasta in a directory by itself
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

        # lineage workflow
        #if not file_exists(op.join(outdir, "completeness.tsv")):
            # clear any existing output directory and create
            # shutil.rmtree(outdir, ignore_errors=True)
            # safe_makedir(outdir)
            # run lineage workflow
        logger.info("Running lineage workflow on %s" % fasta)
        
        completeness_tsv = op.join(outdir, "completeness.tsv")
        
        cmd = "checkm lineage_wf -f {outfile} --tab_table -q -x fasta -t {cores} {binpath} {outdir}".format(outfile=completeness_tsv, outdir=outdir, cores=cores, binpath=bindir)
        print(cmd)
        
        #out = get_stderr(cmd)
        #print(out)
        run(cmd)
        assert op.exists(completeness_tsv)
        #try:
        #    for l in reader(outdir + "/completeness.tsv", header=True, sep="\t"):
        #        completeness = l['Completeness']
        #        break
        #except IOError:
        #    logger.warning("Lineage workflow failed for %s" % fasta)
        #    pass
        completeness = pd.read_csv(completeness_tsv, sep="\t", header=0)
    return completeness
    #return out


def checkm_completeness(sagfile, outfile, cores=10):
    '''calculate checkM completeness value given the SAGs listed in a file

    Args:
        sagfile (str): path to file containing a list of paths, one per line, of SAG fasta files to analyze
        outfile (str): path to location where output table will be written
    Returns:
        tab-delimited file of checkm completeness per SAG
    '''
    logger.info("gathering checkM completeness values for all SAGs listed in file: saglist")
    df = pd.DataFrame(columns=['Bin Id', 'Marker lineage', '# genomes', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity', 'total_bp'])

    saglist = open(sagfile).read().split("\n")
    if len(saglist) < 1:
        return "Error, no SAGs found in sag file"

    for s in saglist:
        if op.exists == False:
            logger.error("SAG not found for %s" % s)
            continue
        if op.isfile == False:
            logger.error("%s is not a file" % s)
            continue
        #sagname = "_".join(op.basename(s).split(".")[:-1])
        #sagpath = s
        completeness = sag_checkm_completeness(s, cores=cores)
        if completeness is None:
            logger.info("completeness stats for %s not determined" % s)
            print("completeness stats for %s not determined" % s)
            continue

        length = count_fasta_bp(s)
        print("sag %s is %s bp in length" % (s, length))

        completeness['total_bp'] = length
        completeness['calculated_length'] = completeness.total_bp * 100/completeness.Completeness
        
        df = pd.concat([df, completeness])
    
    df.to_csv(outfile, sep="\t")
    return outfile


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


@cli.command('completeness', short_help='create dataframe of completeness values for list of sags')
@click.option('--sagfile', help='path to list of SAG files to run through checkm')
@click.option('--outfile', default=None, help='name and path to output csv')
@click.option('--cores', default=10, help='cores to use')
def completeness(sagfile, outfile, cores):
    if op.exists(op.dirname(op.abspath(outfile))) == False:
        logger.warning("Error, cannot find location for output file, using default output scheme name")
        outfile = None

    if outfile is None:
        outfile = ".".join(op.basename(sagfile).split(".")[:-1])+"_completeness.txt"
        logger.info("output file is: %s" % outfile)
        print("output file is called: %s" % outfile)
    else:
        logger.info("output file is called: %s" % outfile)

    outfile = op.abspath(outfile)
    
    out = checkm_completeness(sagfile, outfile, cores)
    return out


if __name__ == '__main__':
    cli()

    
