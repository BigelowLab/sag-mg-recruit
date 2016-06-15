from __future__ import print_function
from Bio import SeqIO
import os.path as op
import logging
import click

'''Convert rRNA regions of genome to 'N's in output fasta file'''

logger = logging.getLogger(__name__)

def mask_sag(input_gb, out_fasta):
    '''output masked contigs with rRNA changed to 'N's
    Args:

    '''
     if input_gb.endswith(".gb") == False or input_gb.endswith(".gbk") == False:
        logger.error("input file does not appear to be in genbank format.  Please check.")
        return None

    with open(input_gb, "rU") as input_handle, open(out_fasta, "w") as oh:
        rrna_count = 0
        for r in SeqIO.parse(input_handle, "genbank"):
            print(">", r.name, sep="", file=oh)
            s = r.seq
            cloc = []
            masked = ""
            for f in r.features:
                if f.type == "rRNA":
                    #if "16S" in str(f.qualifiers['gene']).upper() or "23S" in str(f.qualifiers['gene']).upper():
                    #    print(f)
                    cloc.append(f.location)    # if the 'type' is rRNA, it should be masked... don't have to check for 16 or 23S
                    logger.info('rRNA gene found on contig %s' % r.name)
                    rrna_count += 1

            # if the contig contains one rRNA gene (most common if rRNA present)
            if len(cloc) == 1:
                masked += s[0:cloc[0].start-1]
                masked += 'N'*(cloc[0].end - cloc[0].start)
                masked += s[cloc[0].end-1:]
            # probably won't be instances where below is true
            elif len(cloc) > 1:
                for i in range(0, len(cloc)):
                    # if it's the first entry
                    if i == 0:
                        masked += s[0:cloc[i].start-1]
                        masked += 'N'*(cloc[i].end - cloc[i].start)
                    # if it's the last entry
                    elif i == len(cloc)-1:
                        masked += s[cloc[i-1].end-1:cloc[i].start-1]
                        masked += 'N'*(cloc[i].end - cloc[i].start)
                        masked += s[cloc[i].end:]
                    else:
                        masked += s[cloc[i-1].end-1:cloc[i].start-1]
                        masked += 'N'*(cloc[i].end - cloc[i].start)
            # if no rRNA on contig, just return the sequence, unmasked
            else:
                masked = s

            for i in range(0, len(masked), 80):
                print(masked[i:i+80], file=oh)
    logger.info('%s rRNA genes found in %s' % (rrna_count, op.basename(input_gb))
    return out_fasta


