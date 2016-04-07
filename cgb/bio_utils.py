"""Miscellaneous bioinformatics utility functions."""

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
