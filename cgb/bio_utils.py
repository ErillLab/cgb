"""Miscellaneous bioinformatics utility functions."""

from subprocess import PIPE, Popen

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def complement(seq):
    """Returns the complement of the seq.

    Args:
        seq (string): the DNA sequence.
    Returns:
        string: the complement sequence
    """
    return str(Seq(seq, generic_dna).complement())


def reverse_complement(seq):
    """Returns the reverse complement of the given sequence.

    Args:
        seq (string): the DNA sequence.
    Returns:
        string: the reverse complement sequence
    """
    return str(Seq(seq, generic_dna).reverse_complement())


def weblogo(seqs, filename):
    """Generates the sequence logo for the given sequences.

    Uses WebLogo program (http://weblogo.threeplusone.com/).
    """
    # Sequences in FASTA format
    fasta = '\n'.join('>seq%d\n%s' % (i, seq) for i, seq in enumerate(seqs))
    p = Popen(['weblogo',
               '--format', 'png',
               '--fout', filename,
               '--color-scheme', 'classic',
               '--errorbars', 'YES'],
              stdout=PIPE, stderr=PIPE, stdin=PIPE, close_fds=True)
    p.communicate(input=fasta)
