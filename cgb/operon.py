class Operon:
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
            loc_end = min(self.genome.length, self.first_gene.end-up)

        return self.genome.subsequence(loc_start, loc_end)

    def regulation_probability(self, motif, prior_m):
        """Returns the probability of regulation of the operon.

        Args:
            motif (BindingModel): the model that is used to score the promoter.
            prior_m (float): the prior probability of regulation of the operon.
        """
        random_seqs = self.genome.random_seqs(length=motif.length, count=10**5)
        estimator = motif.bayesian_estimator(bg_scores=random_seqs,
                                             prior_m=prior_m)
        promoter = self.promoter_region()
        promoter_pssm_scores = [motif.score_seq(promoter[i:i+motif.length])
                                for i in xrange(len(promoter)-motif.length+1)]
        return estimator(promoter_pssm_scores, prior_m)

    def __repr__(self):
        return str([g.locus_tag for g in self._genes])