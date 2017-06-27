import pysam
import glob
import os.path as op
import logging
import click
import sys

from scgc.utils import file_transaction, safe_makedir, run, tmp_dir, pigz_file

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger(__name__)

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

def separate_bam(bam, alignedbam, unalignedbam, pctid=95, minlen=150, overlap=0):
    ''' filters bam file based on pctid, minlen and overlap values
    Args:
        bam (string): input bam file
        alignedbam (string): path to aligned bam
        unalignedbam (string): path to unaligned bam
        pctid (numeric): percent identity of read to count as aligned
        minlen (numeric): minimum alignment length to be considered
        overlap (numeric): percent of read to overlap reference to be counted
    Returns:
        tuple of aligned and unaligned bam outputs

    >> bam = 'Tara_test1_vs_Simons_LoCos_Conc.bam'
    >> alignedbam = 'test_aligned.bam'
    >> unalignedbam = 'test_unaligned.bam'
    >> separate_bam(bam, alignedbam, unalignedbam)
    '''
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as ih, \
        pysam.AlignmentFile(alignedbam, "wb", template=ih) as oh_aligned, \
        pysam.AlignmentFile(unalignedbam, "wb", template=ih) as oh_unaligned:
        good = 0
        good_bp = 0
        total = 0
        for i, l in enumerate(ih):
            if l.is_duplicate:
                continue
            elif read_overlap_pctid(l, pctid, minlen, overlap):
                good += 1
                good_bp += l.query_alignment_length
                oh_aligned.write(l)
            else:
                oh_unaligned.write(l)

    return alignedbam, unalignedbam

def extract_fastq(bam, out_fastq):
    ''' Uses bedtools bamtofastq function to extract reads from bam
    Args:
        bam (string): path to bam alignment file
        out_fastq (string): output fastq to write to
    Returns:
        out_fastq (string): path to written output

    >> bam = 'Tara_test1_vs_Simons_LoCos_Conc.pctid95.overlap0.minlen100.bam'
    >> outfastq = 'testout.fastq'
    >> extract_fastq(bam, outfastq) == out_fastq
    '''
    with file_transaction(out_fastq) as temp_oh:
        cmd = "bedtools bamtofastq -i {bam} -fq {fastq}".format(bam=bam, fastq=temp_oh)
        run(cmd)
    return out_fastq

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('smr_dir')
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
def main(smr_dir, pctid, minlen, overlap):
    '''
    >> smr_dir = '/mnt/stepanauskas_nfs/julia/fragdev/extract_reads/'
    >> bam_files(smr_dir) == ['/mnt/stepanauskas_nfs/julia/fragdev/extract_reads/coverage/Tara_test1_vs_Simons_LoCos_Conc.bam']
    >> main(smr_dir)
    '''
    # grab original bam alignment from smr (determined as the ones without .pctidXX suffix)
    bam_files = lambda smr_dir: [i for i in glob.glob(op.join(smr_dir, 'coverage','*.bam')) if 'pctid' not in i]
    for bam in bam_files(smr_dir):
        aligned_out = '{}_aligned.fastq'.format(bam.split(".")[0])
        unaligned_out = '{}_unaligned.fastq'.format(bam.split(".")[0])
        with tmp_dir() as out_dir:
            logger.info("separating aligned and unaligned reads for {}".format(bam))
            abam = op.join(out_dir, "aligned.bam")
            ubam = op.join(out_dir, "unaliged.bam")
            aligned_bam, unaligned_bam = filter_bam(bam, abam, ubam, pctid=pctid, minlen=minlen, overlap=overlap)

            aligned_out = extract_fastq(aligned_bam, aligned_out)
            unaligned_out = extract_fastq(unaligned_bam, unaligned_out)
        logger.info('fastq files created for {}'.format(bam))


if __name__=='__main__':
    main()
