"""The operon module."""

from random import random


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
    def __init__(self, genes, operon_id):
        assert len(genes) > 0
        assert all(gene.strand == genes[0].strand for gene in genes)
        assert all(gene.chromid == genes[0].chromid for gene in genes)
        self._genes = sorted(genes, key=lambda g: g.start)
        self._operon_id = operon_id

    @property
    def chromid(self):
        """Returns the chromosome/plasmid that the operon belongs to."""
        return self.first_gene.chromid

    @property
    def genes(self):
        """Returns the list of genes of the operon."""
        return self._genes

    @property
    def operon_id(self):
        """Returns the operon id"""
        return self._operon_id

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

    @property
    def regulation_probability(self):
        return self.first_gene.regulation_probability

    @property
    def is_probably_regulated(self):
        """Discretizes the regulation probability and returns True or False,
        randomly chosen and weighted by the posterior probability of
        regulation.

        For example, if the P(regulation)=0.9, this method returns True (the
        gene is regulated) 9 out of 10 times.
        """
        return random() <= self.regulation_probability

    def __repr__(self):
        return str([g.locus_tag for g in self._genes])
