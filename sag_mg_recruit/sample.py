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


@cli.command('sample', short_help='sample reads using seqtk')
@click.option('--n', type=click.INT, help='number of reads to subsample')
@click.option('--fastq', help='fastq to subsample')
@click.option('--outfile', help="name of output file")
def sample(fastq, outfile, n):
    out = run_seqtk_sample(fastq, outfile, n)
    return out

if __name__ == '__main__':
    cli()