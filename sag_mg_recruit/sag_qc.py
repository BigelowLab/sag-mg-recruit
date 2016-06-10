from __future__ import print_function

import os.path as op
from scgc.utils import tmp_dir, safe_makedir, run, file_exists
import logging
import pandas as pd
import shutil
import click


logger = logging.getLogger(__name__)

'''
sag qc information for read recruitment
'''

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """completeness via checkm."""
    pass


def sag_checkm_completeness(fasta,  cores):
    logger.info("Running checkm on %s " % fasta)

    prefix = op.basename(fasta).split(".")[0]
    fasta = op.abspath(fasta)
    completeness = 0.0
    number_unique_markers = 0
    number_multi_copy = 0
    taxonomy_contained = "NA"
    taxonomy_sister_lineage = "NA"

    # we want the fasta in a directory by itself
    with tmp_dir() as bindir:
        tmp_fasta = op.join(bindir, op.basename(fasta))
        shutil.copy(fasta, tmp_fasta)
        outdir = op.join(bindir, "outdir")

        # lineage workflow
        if not file_exists(op.join(outdir, "completeness.tsv")):
            # clear any existing output directory and create
            shutil.rmtree(outdir, ignore_errors=True)
            safe_makedir(outdir)
            # run lineage workflow
            logger.info("Running lineage workflow on %s" % fasta)
            run(("checkm lineage_wf -f {outdir}/completeness.tsv --tab_table "
                 "-q -x fasta -t {cores} {binpath} {outdir}").format(outdir=outdir,
                 cores=cores, binpath=bindir))
        #try:
        #    for l in reader(outdir + "/completeness.tsv", header=True, sep="\t"):
        #        completeness = l['Completeness']
        #        break
        #except IOError:
        #    logger.warning("Lineage workflow failed for %s" % fasta)
        #    pass
        completeness = pd.read_csv(op.join(outdir, "completeness.tsv"), sep="\t", header=0)
    return completeness


def checkm_completeness(sagfile, outfile, cores=10):
    '''calculate checkM completeness value given the SAGs listed in a file

    Args:
        sagfile (str): path to file containing a list of paths, one per line, of SAG fasta files to analyze
        outfile (str): path to location where output table will be written
    Returns:
        tab-delimited file of checkm completeness per SAG
    '''
    logger.info("gathering checkM completeness values for all SAGs listed in file: saglist")
    df = pd.DataFrame(columns=['Bin Id', 'Marker lineage', '# genomes', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity'])

    saglist = open(sagfile).read().split("\n")
    for s in saglist:
        #sagname = "_".join(op.basename(s).split(".")[:-1])
        #sagpath = s
        completeness = sag_checkm_completeness(s, cores=cores)
        df = pd.concat([df, completeness])
    
    df.to_csv(outfile, sep="\t")
    return outfile


@cli.command('completeness', short_help='create dataframe of completeness values for list of sags')
@click.option('--sagfile', help='path to list of SAG files to run through checkm')
@click.option('--outfile', help='fastq to subsample')
@click.option('--cores', default=10, help='cores to use')
def sample(sagfile, outfile, cores):
    out = checkm_completeness(sagfile, outfile, cores)
    return out


if __name__ == '__main__':
    cli()

    
