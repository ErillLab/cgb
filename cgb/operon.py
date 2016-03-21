"""The operon module."""


class Operon:
    """Definition for Operon class.

    An operon is a functional unit that is under the control of a single
    promoter.

    For analysis of bacterial gene regulation, it is the basic class which
    contains attributes of genes as well as operon start/end position and
    strand. The class also includes a method to extract promoter region for the
    operon.

    Finally, the class has a method to compute the binding probability of the
    TF to the promoter region. See binding_model.py for details on the Bayesian
    estimator the regulation probability.
    """
    def __init__(self, genes):
        assert len(genes) > 0
        assert all(gene.strand == genes[0].strand for gene in genes)
        assert all(gene.chromid == genes[0].chromid for gene in genes)
        self._genes = sorted(genes, key=lambda g: g.start)

    @property
    def chromid(self):
        """Returns the chromosome/plasmid that the operon belongs to."""
        return self.first_gene.chromid

    @property
    def genes(self):
        """Returns the list of genes of the operon."""
        return self._genes

    @property
    def genome(self):
        """Returns the genome that the operon belongs to."""
        return self.chromid.genome

    @property
    def strand(self):
        """Returns the strand of the operon."""
        return self._genes[0].strand

    @property
    def is_forward_strand(self):
        """Returns true if the operon is on the forward strand."""
        return self.strand == 1

    @property
    def start(self):
        """Returns the start position of the operon."""
        return min(gene.start for gene in self._genes)

    @property
    def end(self):
        """Returns the end position of the operon."""
        return max(gene.end for gene in self._genes)

    @property
    def first_gene(self):
        """Returns the first gene in the operon."""
        if self.is_forward_strand:
            return self._genes[0]
        else:
            return self._genes[-1]

    def promoter_region(self, up=-300, down=+50):
        """Returns the promoter region of the operon."""
        if self.is_forward_strand:
            loc_start = max(0, self.first_gene.start+up)
            loc_end = self.first_gene.start+down
        else:
            loc_start = self.first_gene.end-down
            loc_end = min(self.chromid.length, self.first_gene.end-up)

        return self.chromid.subsequence(loc_start, loc_end)

    def calculate_regulation_probability(self, prior_regulation):
        """Returns the probability of regulation of the operon.
        This is the posterior probability of regulation given the string
        of TF-binding model scores along the promoter region mapping to the operon.

        Args:
            prior_regulation (float): the prior probability of regulation
        """
        #get the TF-binding model adapted to the genome to which the operon belongs
        binding_model = self.genome.TF_binding_model
        #invoke the binding_model method that returns the binding probability for
        #promoter regions assigned to this operon, using the provided prior
        self._regulation_probability = binding_model.binding_probability(
            self.promoter_region(), prior_regulation)
        return self._regulation_probability

    @property
    def regulation_probability(self):
        """Returns the regulation probability of the operon."""
        return self._regulation_probability

    def __repr__(self):
        return str([g.locus_tag for g in self._genes])
