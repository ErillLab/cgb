from Bio.motifs.matrix import PositionWeightMatrix

import bio_utils
import misc


class BindingModel:
    """Class definition for TF-binding model.

    A TF-binding model is constructed from the collection of PWMs and their
    associated weights.
    """
    def __init__(self, collections, weights,
                 background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
        self._background = background
        self._pwm = self._combine_pwms([c.pwm for c in collections], weights)
        self._collection_set = collections

    @property
    def pwm(self):
        """Returns the position-weight-matrix."""
        return self._pwm

    @property
    def length(self):
        """Returns the length of the combined PWM."""
        return self._pwm.length

    @property
    def pssm(self):
        """Returns the position-specific scoring matrix."""
        return self._pwm.log_odds(self.background)

    @property
    def background(self):
        """Returns the background distribution of the model."""
        return self._background

    @property
    def IC(self):
        """Returns the information content of the PSSM."""
        return self.pssm.mean()

    @property
    def patser_threshold(self):
        """Returns the threshold as used in Hertz, Stormo 1999.

        Patser-threshold satisfies the equality between the -log of the
        false-positive-rate and the information-content -- -log(FPR) = IC
        """
        dist = self.pssm.distribution(precision=10**3)
        return dist.threshold_patser()

    @property
    def sites(self):
        """Returns the binding sites of the motifs used to build the model."""
        return [site for coll in self._collection_set for site in coll.sites]

    def score_seq(self, seq):
        """Returns the PSSM score a given sequence.

        The scores from both strands are combined with the soft-max function.
        """
        assert self.length == len(seq)
        complement_seq = bio_utils.complement(seq)
        pssm_score = sum(misc.log2(2**self.pssm[seq[i]][i] +
                                   2**self.pssm[complement_seq[i]][i])
                         for i in range(len(seq)))
        return pssm_score

    def self_score(self):
        """Returns the list of scores of the sites that the model has."""
        return [self.score_seq(site) for site in self.sites]

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
