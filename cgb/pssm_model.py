import math

from cached_property import cached_property
from Bio.Seq import Seq
from Bio.motifs.matrix import PositionWeightMatrix


from .binding_model import TFBindingModel
from .misc import log2
from .misc import temp_file_name
from .bio_utils import weblogo


class PSSMModel(TFBindingModel):
    """Class definition for PSSM model for TF-binding analysis.

    PSSMModel class based on the ubiquitous PSSM model for TF-binding
    analysis. The PSSM method assumes positional independence in the TF-binding
    motif and computes a sum-of-log-likehood ratios as a reasonable
    approximation to the TF-binding energy contribution of any given
    sequence. The likelihood ratio is based on a position-independent
    probability model (the PSWM) and a background model. To make the model
    generic, a uniform background is assumed by default.

    The PSSMModel subclass incorporates a constructor based on the
    weighted mixing of collections of PSWMs that allows instantiating
    species-specific PSSM models based on available evidence in different
    species and a phylogenetic model of these species relationship with the
    target species.
    """
    def __init__(self, collections, weights,
                 background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}, alphabet="ACGT"):
        """Constructor for the PSSMModel class."""
        super(PSSMModel, self).__init__(collections, background, alphabet)
        self._pwm = self._combine_pwms([c.pwm for c in collections], weights, alphabet)

    @cached_property
    def pwm(self):
        """Returns the position-weight-matrix."""
        return self._pwm

    @cached_property
    def length(self):
        """Returns the length of the combined PWM."""
        return self.pwm.length

    @cached_property
    def pssm(self):
        """Returns the position-specific scoring matrix."""
        return self._pwm.log_odds(self.background)

    @cached_property
    def rev_comp_pssm(self):
        """Returns the reverse complement of the PSSM."""
        return self.pssm.reverse_complement()
    
    @cached_property
    def IC(self):
        """Returns the information content of the PSSM."""
        return self.pssm.mean()

    @cached_property
    def patser_threshold(self):
        """Returns the threshold as used in Hertz, Stormo 1999.

        Patser-threshold satisfies the equality between the -log of the
        false-positive-rate and the information-content -- -log(FPR) = IC
        """
        dist = self.pssm.distribution(precision=10**3)
        return dist.threshold_patser()

    def threshold(self, threshold_type='patser'):
        if threshold_type == 'patser':
            thr = self.patser_threshold
        else:
            raise ValueError
        return thr

    @cached_property
    def sites(self):
        """Returns the binding sites of the motifs used to build the model."""
        return [site for coll in self._collection_set for site in coll.sites]

    def score_self(self):
        """Returns the list of scores of the sites that the model has."""
        return [self.score_seq(site) for site in self.sites]

    def _calculate(self, pssm, sequence):
        """Calculates the PSSM score of the sequence with ambiguous letters.

        By default, Biopython's C module is used for scoring. This method is
        called only for sequences that have ambiguous characters.
        """
        score = 0.0
        for pos in xrange(len(sequence)):
            try:
                score += pssm[sequence[pos]][pos]
            except KeyError:
                # TODO(sefa): Handle ambiguous letters other than 'N'
                score += sum(pssm[l][pos] for l in 'ACGT') / 4
        return score

    def score_seq(self, seq, both=True):
        """Returns the PSSM score for a given sequence for all positions.

        The scores from both strands are combined with the soft-max function.

        Args:
            seq (string): the sequence to be scored
        Returns:
            [float]: list of scores of all positions.
        """
        seq = Seq(seq, self.alphabet)
        scores = self.pssm.calculate(seq)
        rc_scores = self.rev_comp_pssm.calculate(seq)

        if self.length == len(seq):
            # Biopython returns single number if len(seq)==len(pssm)
            scores, rc_scores = [scores], [rc_scores]

        # Biopython doesn't handle ambiguous bases well. Calculate score for
        # sites with ambiguous letters.
        for i in xrange(len(scores)):
            if math.isnan(scores[i]):
                site = seq[i:i+self.length]
                scores[i] = self._calculate(self.pssm, site)
                rc_scores[i] = self._calculate(self.rev_comp_pssm, site)

        if both:
            scores = [log2(2**score + 2**rc_score)
                      for score, rc_score in zip(scores, rc_scores)]

        return scores

    @property
    def weblogo_from_pwm(self):
        s = len(self.sites)
        alpha = self._alphabet
        cols = []
        for i in range(self.length):
            cols.append("")
            for let in alpha:
                cols[i] = cols[i] + let * int(self.pwm[let][i] * s)
            if len(cols[i]) < s:
                # column too short
                cols[i] += (alpha * s)[:(s - len(cols[i]))]

        insts = []
        for i in range(s):
            inst = ""
            for j in range(self.length):
                inst += cols[j][i]
            insts.append(inst)

        filename = temp_file_name()
        weblogo(insts, filename)
        return filename

    @staticmethod
    def _combine_pwms(pwms, weights, alphabet):
        """Combines the given PWMs according to the given weights."""
        len = pwms[0].length
        # Check if all PWMs are of the same length.
        assert all(len == pwm.length for pwm in pwms)
        # Combine all PWMs according to given weights
        pwm_vals = {let: [sum(pwm[let][i]*w for pwm, w in zip(pwms, weights))
                          for i in xrange(len)]
                    for let in alphabet}  
        return PositionWeightMatrix(alphabet, pwm_vals)
