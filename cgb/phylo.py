"""Methods for multiple sequence alignment and tree construction."""
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio import Phylo as BioPhylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from cached_property import cached_property

from misc import unique


class Phylo:
    """Class for phylogeny.

    The class Phylogeny contains the phylogenetic tree initialized with the
    given collection of proteins. For the construction of the tree following
    are used:

    - Clustal Omega command-line tool for multiple sequence alignment
    - Distance matrix to compute the distance between two proteins
      (e.g. BLOSUM62, identity matrix)
    - Tree construction algorithm: UPGMA (Unweighted Pair Group Method with
      Arithmetic Mean) or NJ (Neighbor Joining)

    The class also provides methods for outputting the built phylogenetic tree,
    such as drawing it as a string as well as exporting it to a Newick file.
    """
    def __init__(self, proteins, distance_model='identity',
                 tree_algorithm='nj'):

    """Initializes a Phylo object.

        Args:
            proteins (list): list of Protein objects
            distance_model (string): see DistanceCalculator.protein_models
            tree_algorithm (string): 'nj' or 'upgma'
        """
        self._proteins = unique(proteins, lambda p: p.accession_number)
        self._distance_model = distance_model
        self._tree_algorithm = tree_algorithm

    @property
    def proteins(self):
        """Returns Protein objects."""
        return self._proteins

    @cached_property
    def alignment(self):
        """Returns the multiple sequence alignment."""
        return self._clustalo()

    def proteins_to_fasta_file(self, filename):
        """Writes proteins to a temporary FASTA file."""
        with open(filename, 'w') as f:
            for protein in self.proteins:
                f.write(protein.to_fasta())

    def _clustalo(self):
        """Performs Clustal-Omega multiple sequence alignment.

        Args:
            proteins (list): List of Protein objects to be aligned
        Returns:
            MultipleSeqAlignment: A Bio.Align.MultipleSeqAlignment object.
        """
        infile = '/tmp/input.fasta'
        outfile = '/tmp/output.aln'
        self.proteins_to_fasta_file(infile)
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

    @cached_property
    def tree(self):
        """Returns a phylogenetic tree constructed from the given alignment."""
        calculator = DistanceCalculator(self._distance_model)
        constructor = DistanceTreeConstructor(calculator, self._tree_algorithm)
        tree = constructor.build_tree(self.alignment)
        return tree

    def draw_ascii(self):
        """Draws the tree in ASCII format."""
        BioPhylo.draw_ascii(self.tree)

    def to_newick(self, filename):
        """Writes the tree to the given file in newick format."""
        BioPhylo.write(self.tree, filename, 'newick')
