from __future__ import print_function

from sample import run_seqtk_sample
import os.path as op
from scgc.fastx import readfx
from scgc.utils import file_transaction
import argparse

def sample_and_derep(fastq, n, outdir):
    prefix = op.join(outdir, op.basename(fastq).split(".")[0])
    outfile1 = "{prefix}_sub{n}.fastq.gz".format(**locals()) 
    outfile2 = "{prefix}_sub{n}_derep.fastq.gz"
    with file_transaction([outfile1, outfile2]) as tx:
        outfile = run_seqtk_sample(fastq, tx[0], n)
        readset = set()
        seqcount=0
        uniquecount=0

        with open(tx[1], "w") as oh:
            for name, seq, qual in readfx(outfile):
                seqcount += 1
                if seq not in readset:
                    uniquecount += 1
                    readset.add(seq)
                    print("@%s\n%s\n+\n%s\n" % (name, seq, qual), file=oh)
    print("Total sequences examined: %s" % seqcount)
    print("Total unique sequences: %s" % uniquecount)
    return outfile2

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='sample reads and dereplicate')
    p.add_argument('--fastq')
    p.add_argument('--n')
    p.add_argument('--outdir')

    args = p.parse_args()
    sample_and_derep(args.fastq, args.n, args.outdir)