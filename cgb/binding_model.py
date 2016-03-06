import logging
from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np
import scipy.stats
from cached_property import cached_property


class TFBindingModel():
    """Definition for base class TFBindingModel.

    The TFBindingModel is an abstract class that supports arbitrary models for
    TF-binding estimated from the data, using a variety of ML paradigms. The
    TFBindingModel defines a scoring function for specific binding of a TF to a
    DNA sequence, together with a Bayesian statistical model for assessing the
    posterior probability of TF regulation for specific sequences.

    As such, it incorporates the following main methods:
    - A default constructor with known binding sites/regions and background
      sequences
    - A scoring function to score sequences
    - A method to report the score of known binding sites/regions
    - A method to define a threshold for score reporting
    - A method to set up the Bayesian estimator
    - A method to report posterior probabilities of TF regulation on sequences

    Its primary specification is the PSSMModel. See pssm_model.py
    """

    __metaclass__ = ABCMeta

    def __init__(self, collections,
                 background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
        """Constructor for the TFBindingModel class.

        Combines multiple PSWMs using phylogenetic weighting to derive a
        species-specific PWM.
        """
        self._background = background
        self._collection_set = collections
        self._bayesian_estimator = None

    @cached_property
    def background(self):
        """Returns the background distribution of the model."""
        return self._background

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
        logging.debug("Building Bayesian estimator.")
        logging.debug("prior_m: %.2f" % prior_m)
        prior_g = 1.0 - prior_m  # prior probability of the background

        # Estimate the mean/std of functional site scores.
        pssm_scores = self.score_self()
        mu_m, std_m = np.mean(pssm_scores), np.std(pssm_scores)
        logging.debug("pssm_scores - mean: %.2f, std: %.2f" % (mu_m, std_m))
        # Estimate the mean/std of the background scores.
        mu_g, std_g = np.mean(bg_scores), np.std(bg_scores)
        logging.debug("bg_scores - mean: %.2f, std: %.2f" % (mu_g, std_g))

        # Distributions
        pdf_m = scipy.stats.distributions.norm(mu_m, std_m).pdf
        pdf_g = scipy.stats.distributions.norm(mu_g, std_g).pdf
        # Likelihoods
        lh_g = lambda scores: pdf_g(scores)
        lh_m = lambda scores: alpha * pdf_m(scores) + (1-alpha) * pdf_g(scores)
        lh_ratio = lambda scores: np.exp(
            np.sum(np.log(lh_g(scores)) - np.log(lh_m(scores))))

        fun = lambda scores: 1 / (1 + lh_ratio(scores) * prior_g / prior_m)
        self._bayesian_estimator = fun

    def binding_probability(self, seq):
        """Returns the probability of binding to the given seq."""
        pssm_scores = self.score_seq(seq)
        return self.bayesian_estimator(pssm_scores)

    # All of the following methods should be overridden in the subclass.
    @abstractmethod
    def threshold():
        """Returns the threshold for site scores.

        It is an abstract method and should be overridden by the subclass.
        """
        pass

    @abstractmethod
    def score_seq():
        """Returns the score of the given sequence.

        It is an abstract method and should be overridden by the subclass.
        """
        pass

    @abstractproperty
    def length():
        """Returns the length of the binding sites that the model identifies.

        It is an abstract method and should be overridden by the subclass.
        """
        pass
