from cached_property import cached_property

from my_logger import my_logger
from protein import Protein
from blast import BlastNoHitFoundException


class Gene:
    def __init__(self, chromid, seq_feature, product_feature=None):
        """Initializes Gene instance with the given Biopython SeqFeature."""
        self._seq_feature = seq_feature
        self._product_feature = product_feature
        self._chromid = chromid

    @cached_property
    def start(self):
        """Returns the start position of the gene."""
        return self._seq_feature.location.start.position

    @cached_property
    def end(self):
        """Returns the end position of the gene."""
        return self._seq_feature.location.end.position

    @cached_property
    def strand(self):
        """Returns the DNA strand that the gene is on."""
        return self._seq_feature.strand

    @cached_property
    def sequence(self):
        """Returns the gene sequence."""
        return self.chromid.subsequence(self.start, self.end, self.strand)

    @cached_property
    def length(self):
        """Returns the length of the genome."""
        return self.end - self.start

    @cached_property
    def is_forward_strand(self):
        """Returns true if the gene is on the forward strand."""
        return self.strand == 1

    @property
    def chromid(self):
        """Returns the Chromid that the gene belongs to."""
        return self._chromid

    @property
    def genome(self):
        """Returns the Genome that the gene belongs to."""
        return self.chromid.genome

    @property
    def db_xrefs(self):
        """Returns the list of db_xrefs."""
        return self._seq_feature.qualifiers['db_xref']

    @cached_property
    def name(self):
        """Returns the name of the gene."""
        try:
            return self._seq_feature.qualifiers['gene'][0]
        except KeyError:
            return self.locus_tag

    @cached_property
    def locus_tag(self):
        """Returns the locus tag of the gene."""
        locus_tags = self._seq_feature.qualifiers['locus_tag']
        assert len(locus_tags) == 1
        return locus_tags[0]

    @cached_property
    def product_type(self):
        """Returns the product type of the gene."""
        return self._product_feature.type if self._product_feature else ''

    @cached_property
    def product(self):
        """Returns the product. Returns None if no product."""
        try:
            p = self._product_feature.qualifiers['product'][0]
        except (AttributeError, KeyError):
            p = ''
        return p

    @cached_property
    def is_protein_coding_gene(self):
        """Returns true if the gene is a protein coding gene."""
        return self.product_type == 'CDS'

    @cached_property
    def protein_accession_number(self):
        """Returns the accession number of the protein coded by this gene."""
        # TODO(sefa): throw execption if not protein coding gene
        return self._product_feature.qualifiers['protein_id'][0]

    def to_protein(self):
        """Returns the protein object for the gene."""
        return Protein(self.protein_accession_number, self.name)

    def find_homolog_in_genome(self, genome):
        """Returns the homologous gene in the given genome."""
        return genome.find_gene_homolog(self)

    def reciprocal_blast_hit(self, genome):
        """Returns the reciprocal best hit of the gene in the given genome."""
        try:
            # Find the best hit in the other genome
            best_hit, _ = self.find_homolog_in_genome(genome)
            # Check if it is the best reciprocal hit
            reciprocal_hit, _ = best_hit.find_homolog_in_genome(self.genome)
            if self == reciprocal_hit:
                return best_hit
        except BlastNoHitFoundException:
            my_logger.debug("No reciprocal BLAST hit for %s." % self.locus_tag)
            return None

    def distance(self, other):
        """Returns the distance between two genes.

        Used to compute the intergenic distance between two neighboring genes.
        """
        return max(self.start, other.start) - min(self.end, other.end)

    def to_fasta(self):
        """Returns the gene sequence as a string in FASTA format."""
        description = '>%s\n' % self.locus_tag
        line_length = 60
        sequence = '\n'.join(self.sequence[i:i+line_length]
                             for i in range(0, self.length, line_length))
        return description + sequence

    def __repr__(self):
        return self.locus_tag
