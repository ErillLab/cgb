"""Methods for multiple sequence alignment and tree construction."""
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO


def proteins_to_fasta_file(proteins, filename):
    """Writes proteins to a temporary FASTA file."""
    with open(filename, 'w') as f:
        for protein in proteins:
            f.write(protein.to_fasta())


def clustalo(proteins, infile='/tmp/input.fasta', outfile='/tmp/output.aln'):
    """Performs Clustal-Omega multiple sequence alignment.

    Args:
        proteins (list): List of Protein objects to be aligned
        infile (string): Clustal input file name
        outfile (string): Clustal output file name

    Returns:
        MultipleSeqAlignment: A Bio.Align.MultipleSeqAlignment object.
    """
    proteins_to_fasta_file(proteins, infile)
    clustalo_cline = ClustalOmegaCommandline(
        'clustalo',             # executable
        infile=infile,          # input file name
        outfile=outfile,        # output file name
        outfmt='clustal',       # output format
        verbose=True,           # verbose output
        auto=True,              # set options automatically
        force=True)             # force file overwriting
    stdout, stderr = clustalo_cline()
    print stderr

    align = AlignIO.read(outfile, 'clustal')
    return align
