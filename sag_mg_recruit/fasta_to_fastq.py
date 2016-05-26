#!usr/bin/env python

from scgc.utils import file_transaction, pigz_file
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
import argparse

'''
Create fastq from 454 fasta and qual files
'''

def combine_fasta_qual(fas, qual, outfile, cores=8):
    if outfile.endswith(gz) == False:
        outfile = outfile+".gz"

    with file_transaction(outfile) as tx_out:
        with open(fas) as fin, open(qual) as qin, open(tx_out, "w") as oh:
            for rec in PairedFastaQualIterator(fin, qin):
                SeqIO.write(rec, oh, "fastq")
    outfile = pigz_outfile(outfile, cores)  
    return outfile  


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('fasta')
    p.add_argument('qual')
    p.add_argument('outfile')
    p.add_argument('cores', default=8)
    
    args = p.parse_args()

    combine_fasta_qual(args.fasta, args.qual, args.outfile)

