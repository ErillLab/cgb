from protein import Protein

class Gene:
    def __init__(self, chromid, seq_feature, product_feature=None):
        """Initializes Gene instance with the given Biopython SeqFeature."""
        self._seq_feature = seq_feature
        self._product_feature = product_feature
        self._chromid = chromid

    @property
    def start(self):
        """Returns the start position of the gene."""
        return self._seq_feature.location.start.position

    @property
    def end(self):
        """Returns the end position of the gene."""
        return self._seq_feature.location.end.position

    @property
    def strand(self):
        """Returns the DNA strand that the gene is on."""
        return self._seq_feature.strand

    @property
    def sequence(self):
        """Returns the gene sequence."""
        return self.chromid.subsequence(self.start, self.end, self.strand)

    @property
    def length(self):
        """Returns the length of the genome."""
        return self.end - self.start

    @property
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

    @property
    def name(self):
        """Returns the name of the gene."""
        try:
            return self._seq_feature.qualifiers['gene'][0]
        except KeyError:
            return self.locus_tag

    @property
    def locus_tag(self):
        """Returns the locus tag of the gene."""
        locus_tags = self._seq_feature.qualifiers['locus_tag']
        assert len(locus_tags) == 1
        return locus_tags[0]

    @property
    def product_type(self):
        """Returns the product type of the gene."""
        return self._product_feature.type if self._product_feature else ''

    @property
    def product(self):
        """Returns the product. Returns None if no product."""
        try:
            p = self._product_feature.qualifiers['product'][0]
        except (AttributeError, KeyError):
            p = ''
        return p

    @property
    def is_protein_coding_gene(self):
        """Returns true if the gene is a protein coding gene."""
        return self.product_type == 'CDS'

    @property
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
