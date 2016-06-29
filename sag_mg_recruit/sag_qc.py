from __future__ import print_function

import os.path as op
from scgc.utils import tmp_dir, safe_makedir, run, file_exists
import logging
import pandas as pd
import shutil
import click
from itertools import groupby
import gzip
from Bio import SeqIO


logger = logging.getLogger(__name__)

'''
sag transformations and qc information for read recruitment
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
    if op.isdir == True or op.exists == False:
        return None
    
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

        logger.info("Running lineage workflow on %s" % fasta)
        
        completeness_tsv = op.join(outdir, "completeness.tsv")
        
        cmd = "checkm lineage_wf -f {outfile} --tab_table -q -x fasta -t {cores} {binpath} {outdir}".format(outfile=completeness_tsv, outdir=outdir, cores=cores, binpath=bindir)
        
        logger.info("running checkm lineage, command is: {cmd}".format(**locals()))
        run(cmd)
        completeness = pd.read_csv(completeness_tsv, sep="\t", header=0)
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
        completeness = sag_checkm_completeness(s, cores=cores)
        if completeness is None:
            logger.info("completeness stats for %s not determined" % s)
            #print("completeness stats for %s not determined" % s)
            continue

        length = count_fasta_bp(s)
        #print("sag %s is %s bp in length" % (s, length))

        completeness['total_bp'] = length
        completeness['calculated_length'] = int(completeness.total_bp * 100/completeness.Completeness)
        
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


def mask_sag(input_gb, out_fasta):
    '''output masked contigs with rRNA changed to 'N's
    Args:
        input_gb (str): path to annotated input genbank formatted genome
        out_fasta (str): where to write the output fasta to
    Returns:
        fasta file with rRNA regions masked with 'N's
        out_fasta (str)
    '''
    #if input_gb.endswith(".gb") == False or input_gb.endswith(".gbk") == False:
    #    logger.error("input file does not appear to be in genbank format.  Please check.")
    #    return None

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
                        "23S" in str(f.qualifiers['gene']).upper())):
                            cloc.append(f.location)    # if the 'type' is rRNA, it should be masked... don't have to check for 16 or 23S
                            logger.info('rRNA gene found on contig %s' % r.name)
                            rrna_count += 1      
                    elif ('product' in f.qualifiers and 
                         ("16S" in str(f.qualifiers['product']).upper() or 
                        "23S" in str(f.qualifiers['product']).upper())):
                        #print(f)
                            cloc.append(f.location)    # if the 'type' is rRNA, it should be masked... don't have to check for 16 or 23S
                            logger.info('rRNA gene found on contig %s' % r.name)
                            rrna_count += 1
                    else:
                        print(f)

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
    logger.info('%s rRNA genes found in %s' % (rrna_count, op.basename(input_gb)))
    return out_fasta


@cli.command('completeness', short_help='create dataframe of completeness values for list of sags')
@click.option('--sagfile', help='path to list of SAG files to run through checkm')
@click.option('--outfile', default=None, help='name and path to output csv')
@click.option('--cores', default=10, help='cores to use')
@click.option('--log', default="completeness_log.txt", help='name of logfile')
def completeness(sagfile, outfile, cores, log):
    logging.basicConfig(filename=log, level=logging.DEBUG)
    
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


@cli.command('mask', short_help='mask rRNA genes in sag genome')
@click.option('--input_gb', help='path to input annotated genbank file')
@click.option('--out_fasta', default=None, help='where to write masked output fasta file to')
@click.option('--log', default=None, help='name of logfile')
def run_mask_sag(input_gb, out_fasta, log):
    if log == None:
        log = "mask_%s.log" % "_".join(input_gb.split(".")[:-1])
        print("Logfile is: %s" % log)
    logging.basicConfig(filename=log, level=logging.DEBUG)
    if out_fasta == None:
        out_fasta = op.basename(input_gb).split(".")[0]+"_masked.fasta"
        logger.info("output fasta file will be: %s" % out_fasta)
        print("output fasta file will be: %s" % out_fasta)
    out = mask_sag(input_gb, out_fasta)    
    return out


if __name__ == '__main__':
    cli()
