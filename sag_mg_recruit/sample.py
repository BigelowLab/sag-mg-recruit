from __future__ import print_function
import click
import logging

from scgc.utils import file_exists, file_transaction, run, pigz_file

logger = logging.getLogger(__name__)

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """read recruitment."""
    pass


@cli.command('sample', short_help='sample reads using seqtk')
@click.option('--n', type=click.INT, help='number of reads to subsample')
@click.option('--fastq', help='fastq to subsample')
@click.option('--outfile', help="name of output file")
def run_seqtk_sample(fastq, outfile, n, seed=37):
    """Subsample incoming paired-end fastqs to `n` reads (serially).

    Args:
        fastqs (str): path to fastq
        outfile (str): path of output fastq paths; output files are always gzipped
        n (int): number of subsampled reads
        seed (int): for random selection of reads

    Returns:
        str: subsampled reads file path
    """
    if file_exists(outfile):
        return outfile

    logger.info("Subsampling to %d reads" % n)
    with file_transaction(outfile) as tx:
        cmd = "seqtk sample -s {seed} {fastq} {number} | gzip > {out}".format(seed=seed, fastq=fastq, number=n, out=tx)
        run(cmd)
    print("%s created" % outfile)
    return outfile


def run_kmernorm(fastq, out_file, options, cores=1):
    """
    perform digital normalization using kmernorm. unzip input file if gzipped
    and gzip input file after kmernorm.

    parameters
        fastq : file path; gzipped input will be decompressed and recompressed
        out_file : file path; file ending with '.gz' will be compressed
        options : kmernorm options
        cores : cpu cores to utilize

    returns
        normalized file path and surviving read pairs : tuple (string, int)
    """
    if file_exists(out_file):
        return out_file

    gzip_result = True if out_file.endswith('.gz') else False
    out_file = out_file.rsplit('.gz', 1)[0]
    if fastq.endswith(".gz"):
        fastq = gunzip_file(fastq)
    with file_transaction(out_file) as tx:
        run("kmernorm %s %s > %s" % (options, fastq, tx),
            description="Normalizing with kmernorm")
    if gzip_result:
        out_file = pigz_file(out_file, cores)
    # compress input file regardless
    fastq = pigz_file(fastq, cores)
    return out_file


if __name__ == '__main__':
    cli()