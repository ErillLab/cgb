"""Miscellaneous bioinformatics utility functions."""


def complement(seq):
    """Returns the complement of the seq.

    Args:
        seq (string): the DNA sequence.
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement_dict[b] for b in seq])


def reverse_complement(seq):
    """Returns the reverse complement of the given sequence.

    Args:
        seq (string): the DNA sequence.
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement_dict[b] for b in reversed(seq)])
