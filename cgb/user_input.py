import json


class UserInput:
    """Class definition for UserInput.

    UserInput class gives an interface to access the run parameters provided by
    the user. It allows us to make changes in the input format (currently,
    JSON) relatively easily. It also checks the content of the input and sets
    default values for parameters not specified by the user.
    """
    def __init__(self, input_filename):
        """Constructs Input object using the input file."""
        self._input = {}
        with open(input_filename) as f:
            for k, v in json.load(f).items():
                self._input[k] = v

    @property
    def genome_name_and_accessions(self):
        """Returns the list of (genome name, accession_numbers) tuples."""
        return [(g['name'], g['accession_numbers'])
                for g in self._input['genomes']]

    @property
    def genome_names(self):
        """Returns the list of genome names used."""
        return [g['name'] for g in self._input['genomes']]

    @property
    def protein_accessions(self):
        """Returns the protein accession numbers."""
        return [p['protein_accession'] for p in self._input['motifs']]

    @property
    def TF_name(self):
        """Returns the TF of interest."""
        return self._input['TF']

    @property
    def sites_list(self):
        """Returns the lists of binding sites."""
        return [m['sites'] for m in self._input['motifs']]

    @property
    def protein_accessions_and_sites(self):
        """Zips protein accessions and binding sites."""
        return zip(self.protein_accessions, self.sites_list)

    @property
    def prior_regulation_probability(self):
        """Returns the prior probability of regulation.

        It is used by the Bayesian estimator of operon regulation as the prior
        probability of an operon being regulated.

        If not provided by the user, it is estimated by predicting operons in
        the genome of the TF for which the motif is provided and dividing the
        number of sites in the motif (assuming that there is one site per
        operon) by the number of operons. For instance, if 30 sites are
        available for LexA in E. coli, then the prior for regulation is
        30/~2300 operons. See "get_prior" in "main.py" for implementation.
        """
        try:
            value = self._input['prior_regulation_probability']
        except KeyError:
            value = None
        return value

    @property
    def has_prior_probability_set(self):
        """Returns true if the prior probability of regulation is provided."""
        return 'prior_regulation_probability' in self._input

    @property
    def probability_threshold(self):
        """Returns the threshold for regulation probabilities.

        Only the operons with a regulation probability above the threshold will
        be reported.

        The default value for the probability threshold is 0.5
        """
        try:
            value = self._input['posterior_probability_threshold_for_reporting']
        except KeyError:
            value = 0.5
        return value

    @property
    def phylogenetic_weighting(self):
        """Returns if the phylogenetic-weighting option is on."""
        try:
            value = self._input['phylogenetic_weighting']
        except KeyError:
            value = False
        return value

    @property
    def site_count_weighting(self):
        """Returns if the site count-weighting option is on."""
        try:
            value = self._input['site_count_weighting']
        except KeyError:
            value = False
        return value

    @property
    def operon_prediction_probability_threshold(self):
        """Returns the threshold for posterior probability of regulation used
        for operon prediction. The default value is 0.5.
        """
        try:
            value = self._input['posterior_probability_threshold_for_operon_prediction']
        except KeyError:
            value = 0.5
        return value
