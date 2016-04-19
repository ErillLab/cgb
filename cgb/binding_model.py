import logging

# Enforces the definition of virtual methods in derived classes
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

    @cached_property
    def background(self):
        """Returns the background distribution of the model."""
        return self._background

    @property
    def site_collections(self):
        """Returns the site collections used to build the model."""
        return self._collection_set

    @property
    def bayesian_estimator(self):
        """Returns the Bayesian estimator for computing P(regulation|data).

        The Bayesian estimator is set at "build_bayesian_estimator" method. It
        is a function that takes a list of scores as the argument, and returns
        the probability of binding.
        """
        return self._bayesian_estimator

    def build_bayesian_estimator(self, bg_scores):
        """Builds a Bayesian estimator of the probability of regulation.

        It computes probability distributions for PSSM scores of binding sites
        and background sequences. The probability distributions are stored for
        later use to compute posterior probability of binding on a given
        sequence.

        Args:
            bg_scores (list): List of PSSM scores of sequences from background.
        """
        logging.debug("Building Bayesian estimator.")
        # Estimate the mean/std of functional site scores.
        pssm_scores = self.score_self()
        self._mu_m, self._std_m = np.mean(pssm_scores), np.std(pssm_scores)
        # Estimate the mean/std of the background scores.
        self._mu_bg, self._std_bg = np.mean(bg_scores), np.std(bg_scores)

    def binding_probability(self, seq, p_motif, alpha=1/350.0):
        """Returns the probability of binding to the given seq.

        Args:
            seq (string): the sequence to compute the probability of TF-binding
            p_motif (float): prior probability of binding
            alpha (float): The mixing ratio.
        Returns:
            float: the probability of TF-binding to the sequence.
        """
        pssm_scores = self.score_seq(seq)
        p_bg = 1.0 - p_motif
        # Probability density functions
        pdf_m = scipy.stats.distributions.norm(self._mu_m, self._std_m).pdf
        pdf_bg = scipy.stats.distributions.norm(self._mu_bg, self._std_bg).pdf
        # Compute the likelihood of the seq from motif and background
        lh_m = (alpha * pdf_m(pssm_scores) + (1-alpha) * pdf_bg(pssm_scores))
        lh_bg = pdf_bg(pssm_scores)
        # Compute the likelihood ratio
        lh_ratio = np.exp(np.sum(np.log(lh_bg) - np.log(lh_m)))
        return 1 / (1 + lh_ratio * p_bg / p_motif)

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
