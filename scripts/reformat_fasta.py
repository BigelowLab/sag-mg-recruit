from __future__ import print_function
from itertools import groupby
import click


def read_fasta(file_handle):
    '''Fasta iterator'''
    for header, group in groupby(file_handle, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('fasta')
@click.argument('new_fasta')
@click.option('--line-len', help='specify the length of the wraped sequence', default=60, show_default=True)
def main(fasta, new_fasta, line_len):
    ''' Rewrites fasta, removing any spaces and new line characters from previous fasta
    Args:
        fasta (path): fasta file you'd like to rewrite
        new_fasta (path): where you'd like to write the reformatted fasta to
    Returns:
        new fasta file
    '''
    with open(new_fasta, "w") as oh, open(fasta) as ih:
        for name, seq in read_fasta(ih):
            print(">{name}".format(name=name), file=oh)

            seq_fix = seq.replace(" ","")
            for i in range(0, len(seq_fix), line_len):
                print(seq_fix[i:i+line_len], file=oh)

if __name__=='__main__':
    main()
