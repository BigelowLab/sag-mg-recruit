#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function
from __future__ import division

import numpy as np
import math
from scgc.utils import file_transaction, multiprocess, pigz_file, file_exists
from scgc.fastx import readfx
import pandas as pd
import multiprocessing
import logging
import click
import numpy as np


logger = logging.getLogger(__name__)

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """read recruitment."""
    pass


class OnlineVariance(object):
    """
    Welford's algorithm computes the sample variance incrementally.
    from: http://stackoverflow.com/questions/5543651/computing-standard-deviation-in-a-stream
    
    Example use:
    > N = 100
    > data = np.random.random(N)
    > ov = OnlineVariance(ddof=0)
    > for d in data:
    >     ov.include(d)
    > std = ov.std
    > print(std)
    """

    def __init__(self, iterable=None, ddof=1):
        self.ddof, self.n, self.mean, self.M2 = ddof, 0, 0.0, 0.0
        if iterable is not None:
            for datum in iterable:
                self.include(datum)

    def include(self, datum):
        self.n += 1
        self.delta = datum - self.mean
        self.mean += self.delta / self.n
        self.M2 += self.delta * (datum - self.mean)
        self.variance = self.M2 / (self.n - self.ddof)

    @property
    def std(self):
        return np.sqrt(self.variance)


def read_fq_string(fq_list):
    for n in range(0, len(fq_list)-1, 4):
        name = fq_list[n].replace("@","")
        seq = fq_list[n+1]
        qual = fq_list[n+3]
        yield name, seq, qual


def chop_read((name, seq, qual), length, minimum):
    #pieces = round(len(seq)/length)
    #div = int(round(len(seq)/pieces))
    #div = int(math.floor(len(seq)/length))
    div = int(round(np.random.normal(loc=length, scale=10)))
    #div = length
    newrecords= ""
    #for i, n in enumerate(range(0, len(seq), div)):
    for i, n in enumerate(range(0, len(seq), div)):
        newname = "%s_%s" % (name, i+1)
        newseq = seq[n:n+div]
        newqual = qual[n:n+div]
        
        #if abs(len(seq[n:len(seq)])-length) <= 30:
        #    newseq = seq[n:len(seq)]
        #    newqual = qual[n:len(seq)]
        #    newrecords +=  "@%s\n%s\n+\n%s\n" % (newname, newseq, newqual)
        #    break

        if len(newseq) > minimum:
            newrecords += "@%s\n%s\n+\n%s\n" % (newname, newseq, newqual)

    return str(newrecords)


def chop_reads(fastq, out_file, length, minlength, cores=10):
    
    if not out_file.endswith('.gz'):
        out_file = out_file + '.gz'
    
    if file_exists(out_file):
        print("output file exists, chop reads will not repeat itself.")
        return out_file
    
    pd.set_option('display.float_format', lambda x: '%.0f' % x)
    p = multiprocessing.Pool(cores)
    
    original_count = 0
    readcount = 0
    bpcount = 0
    size_var = OnlineVariance(ddof=0)
    
    with file_transaction(out_file) as tx_outfile:
        with open(tx_outfile, "w") as txoh:
            for result in multiprocess(chop_read, readfx(fastq), length, minlength, pool=p):
    
                if type(result) == list:
                    for r in result:
                        print(r, file=txoh)
                        original_count += 1
                        lines = r.split("\n")
                        for name, seq, qual in read_fq_string(lines):
                            readcount += 1
                            bpcount += len(seq)
                            size_var.include(len(seq))
                else:
                    print(result, file=txoh)
                    original_count += 1
                    lines = result.split("\n")
                    for name, seq, qual in read_fq_string(lines):

                        readcount += 1
                        bpcount += len(seq)
                        size_var.include(len(seq))
    
    meanbp = bpcount/readcount
    readbp_std = size_var.std
    out_file = pigz_file(out_file, cores)
    
    outdata = pd.Series(data=[original_count, readcount, bpcount, meanbp, readbp_std], 
                        index=['original_count', 'read_count', 'bp_count', 'mean_read_len', 'read_len_std'])
    print(outdata)
    outdata.to_csv(out_file.replace("fastq.gz","chop_data"), sep="\t")
    
    return out_file

@cli.command('chop', short_help='cut reads into smaller reads using ')
@click.option('--length', type=click.INT, help='target length')
@click.option('--fastq', help='fastq to subsample')
@click.option('--minlength', type=click.INT, help="minimum length necessary")
@click.option('--outfile', help="name of output file")
@click.option('--cores', default=10, help="number of cores to use")
def run_chop_reads(length, fastq, minlength, outfile, cores):
    out = chop_reads(fastq, outfile, length, minlength, cores)
    return out

if __name__ == '__main__':
    cli()