from Bio.motifs.matrix import PositionWeightMatrix
import numpy as np
import scipy.stats

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
        self._bayesian_estimator = None

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

    def score_self(self):
        """Returns the list of scores of the sites that the model has."""
        return [self.score_seq(site) for site in self.sites]

    @property
    def bayesian_estimator(self):
        """Returns the Bayesian estimator for computing P(regulation)."""
        return self._bayesian_estimator

    def build_bayesian_estimator(self, bg_scores, prior_m, alpha=1.0/350):
        """Builds a Bayesian estimator of the probability of regulation.

        The attribute '_bayesian_estimator' is set to a function that takes a
        list of PSSM scores of a sequence and returns the probability of
        binding to that sequence.

        Args:
            bg_scores (list): List of PSSM scores of sequences from background.
            prior_m (float): The prior probability of regulation.
            alpha (float): The mixing ratio.
        """
        prior_g = 1.0 - prior_m  # prior probability of the background

        # Estimate the mean/std of functional site scores.
        pssm_scores = self.score_self()
        mu_m, std_m = np.mean(pssm_scores), np.std(pssm_scores)
        # Estimate the mean/std of the background scores.
        mu_g, std_g = np.mean(bg_scores), np.std(bg_scores)
        # Distributions
        pdf_m = scipy.stats.distributions.norm(mu_m, std_m).pdf
        pdf_g = scipy.stats.distributions.norm(mu_g, std_g).pdf
        # Likelihoods
        lh_g = lambda scores: pdf_g(scores)
        lh_m = lambda scores: alpha * pdf_m(scores) + (1-alpha) * pdf_g(scores)
        lh_ratio = lambda scores: np.exp(
            np.sum(np.log(lh_g(scores)) - np.log(lh_m(scores))))

        fun = lambda scores: 1 / (1 + lh_ratio(scores) * prior_m / prior_g)
        self._bayesian_estimator = fun

    def binding_probability(self, seq):
        """Returns the probability of binding to the given seq."""
        pssm_scores = [self.score_seq(seq[i:i+self.length])
                       for i in xrange(len(seq)-self.length+1)]
        return self.bayesian_estimator(pssm_scores)

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
