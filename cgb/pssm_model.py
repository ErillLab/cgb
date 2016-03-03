from cached_property import cached_property
from Bio.Seq import Seq
from Bio.motifs.matrix import PositionWeightMatrix


from binding_model import TFBindingModel
from misc import log2


class PSSMModel(TFBindingModel):
    """Class definition for PSSM model for TF-binding analysis."""
    def __init__(self, collections, weights,
                 background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
        """Constructor for the PSSMModel class."""
        super(PSSMModel, self).__init__(collections, weights, background)
        self._pwm = self._combine_pwms([c.pwm for c in collections], weights)

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
    def reverse_complement_pssm(self):
        """Returns the reverse complement of the PSSM."""
        return self.pssm.reverse_complement()

    @property
    def alphabet(self):
        """Returns the alphabet of the motif."""
        return self.pwm.alphabet

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

    def score_seq(self, seq):
        """Returns the PSSM score for a given sequence for all positions.

        The scores from both strands are combined with the soft-max function.

        Args:
            seq (string): the sequence to be scored
        Returns:
            [float]: list of scores of all positions.
        """
        seq = Seq(seq, self.alphabet)
        scores = self.pssm.calculate(seq)
        rc_scores = self.reverse_complement_pssm.calculate(seq)

        if self.length == len(seq):
            # Biopython returns single number if len(seq)==len(pssm)
            scores, rc_scores = [scores], [rc_scores]

        return [log2(2**score + 2**rc_score)
                for score, rc_score in zip(scores, rc_scores)]

    @staticmethod
    def _combine_pwms(pwms, weights):
        """Combines the given PWMs according to the given weights."""
        len = pwms[0].length
        alphabet = pwms[0].alphabet
        # Check if all PWMs are of the same length.
        assert all(len == pwm.length for pwm in pwms)
        # Check if all PWMs have the same alphabet -- 'ACGT'
        assert all(alphabet == pwm.alphabet for pwm in pwms)
        # Normalize weights
        weights = [float(weight)/sum(weights) for weight in weights]
        # Combine all PWMs according to given weights
        pwm_vals = {let: [sum(pwm[let][i]*w for pwm, w in zip(pwms, weights))
                          for i in xrange(len)]
                    for let in alphabet.letters}
        return PositionWeightMatrix(alphabet, pwm_vals)
