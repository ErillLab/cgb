from cached_property import cached_property

from .my_logger import my_logger
from .protein import Protein
from .blast import BlastNoHitFoundException


class Gene:
    """Definition for Gene class.

    Gene class represents the physical and functional units of the genome:
    genes. The class provides access to primary gene features such as the start
    and end position, the strand that it is lying on.
    """

    def __init__(self, index, chromid, seq_feature, product_feature=None):
        """Initializes Gene instance with the given Biopython SeqFeature."""
        self._index = index
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

    @cached_property
    def operon(self):
        """Return the operon that the gene belongs to."""
        opr, = [o for o in self.chromid.operons if self in o.genes]
        return opr

    @cached_property
    def upstream_gene(self):
        """Returns the gene upstream.

        Returns None if there is no such gene.
        """
        genes = self.chromid.genes
        if self.is_forward_strand and self._index > 0:
            return genes[self._index-1]
        elif not self.is_forward_strand and self._index < len(genes)-1:
            return genes[self._index+1]
        return None

    def upstream_noncoding_region_location(self, up=None, down=50):
        """Returns the upstream non-coding region.

        Args:
            up (int): the relative start position of the upstream region. If
                none, the returned region is up to the previous gene.
            down (int): the relative end position of the returned region.
        Returns: (int, int) start and end position of the region
        """
        def upstream_forward_strand_gene():
            """Returns the upstream region of a gene in forward strand."""
            if up:
                loc_start = max(0, self.start-up)
            elif self.upstream_gene:
                loc_start = min(self.upstream_gene.end-50, self.start)
            else:
                loc_start = 0
            loc_end = self.start + down
            return loc_start, loc_end

        def upstream_reverse_strand_gene():
            """Returns the upstream region of a gene in reverse strand."""
            if up:
                loc_end = min(self.chromid.length, self.end+up)
            elif self.upstream_gene:
                loc_end = max(self.upstream_gene.start+50, self.end)
            else:
                loc_end = self.chromid.length
            loc_start = self.end - down
            return loc_start, loc_end

        return (upstream_forward_strand_gene() if self.is_forward_strand
                else upstream_reverse_strand_gene())

    def upstream_noncoding_region_sequence(self, up=None, down=50):
        """Returns the sequence of the upstream region of the gene."""
        start, end = self.upstream_noncoding_region_location(up, down)
        return self.chromid.subsequence(start, end)

    def promoter_region(self, up=-300, down=+50):
        """Returns the promoter region of the gene."""
        if self.is_forward_strand:
            loc_start = max(0, self.start+up)
            loc_end = self.start+down
        else:
            loc_start = self.end-down
            loc_end = min(self.chromid.length, self.end-up)

        return self.chromid.subsequence(loc_start, loc_end)

    def calculate_regulation_probability(self, prior_regulation):
        """Returns the posterior probability of regulation given the TF-binding
        model scores along the gene promoter region.

        Args:
            prior_regulation (float): the prior probability of regulation
        """
        # Get the TF-binding model adapted to the genome to which the operon
        # belongs.
        binding_model = self.genome.TF_binding_model
        # Invoke the binding_model method that returns the binding probability
        # for promoter regions assigned to this operon, using the provided
        # prior
        self._regulation_probability = binding_model.binding_probability(
            self.promoter_region(), prior_regulation)
        return self._regulation_probability

    @property
    def regulation_probability(self):
        """Returns the regulation probability of the operon."""
        return self._regulation_probability

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
        """Returns the name of the gene.

        If the gene has no name, the locus tag is returned.
        """
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
        """Returns the product type of the gene.

        If the gene has no product information, empty string is returned.
        """
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
        try:
            return self._product_feature.qualifiers['protein_id'][0]
        except KeyError:
            raise NotProteinCodingGeneException

    def to_protein(self):
        """Returns the protein object for the gene."""
        return Protein(self.protein_accession_number)

    def find_homolog_in_genome(self, genome):
        """Returns the homologous gene in the given genome."""
        return genome.find_gene_homolog(self)

    def reciprocal_blast_hit(self, genome, cache):
        """Returns the reciprocal best hit against a provided genome.

        The function uses "find_homolog_in_genome", which calls the
        "find_gene_homolog" in the respective genome to search using BLAST.
        """
        # Check if the cache contains the reciprocal BLAST hit
        key = (self.locus_tag,
               tuple([c.accession_number for c in genome.chromids]))
        if key in cache.keys():
            my_logger.debug("Key found in BLAST cache: %s" % str(key))
        if key not in cache.keys():
            my_logger.debug("Key NOT found in BLAST cache: %s" % str(key))
            # Search for reciprocal best BLAST hit
            cache[key] = None
            try:
                # Find the best hit in the other genome
                best_hit, _ = self.find_homolog_in_genome(genome)
                # Check if it is the best reciprocal hit
                reciprocal_hit, _ = best_hit.find_homolog_in_genome(self.genome)
                # If the reciprocal hit is the very own gene, return the target
                # genome hit as the best_reciprocal BLAST hit for the gene
                if self == reciprocal_hit:
                    cache[key] = best_hit.locus_tag
                    other_key = (best_hit.locus_tag,
                                 tuple([c.accession_number
                                        for c in best_hit.genome.chromids]))
                    cache[other_key] = self.locus_tag
            # Handle lack of significant results from the BLAST search (in
            # either direction)
            except BlastNoHitFoundException:
                my_logger.debug("No reciprocal BLAST hit for %s." % self.locus_tag)
        return genome.get_gene_by_locus_tag(cache[key]) if cache[key] else None

    def distance(self, other):
        """Returns the distance between two genes.

        Used to compute the intergenic distance between two neighboring genes.
        """
        return max(self.start, other.start) - min(self.end, other.end)

    def distance_to_region(self, region_start, region_end):
        """Returns the distance to the region given its location."""
        return max(self.start, region_start) - min(self.end, region_end)

    def to_fasta(self):
        """Returns the gene sequence as a string in FASTA format."""
        description = '>%s\n' % self.locus_tag
        line_length = 60
        sequence = '\n'.join(self.sequence[i:i+line_length]
                             for i in range(0, self.length, line_length))
        return description + sequence

    def __repr__(self):
        return self.locus_tag


class NotProteinCodingGeneException:
    pass
