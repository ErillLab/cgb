from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet.IUPAC import unambiguous_dna


class SiteCollection:
    """The class definition for the collection of binding sites"""
    def __init__(self, sites, TF, pseudocounts=1):
        self._TF = TF
        instances = [Seq(site, unambiguous_dna) for site in sites]
        self._motif = motifs.create(instances)
        self._motif.pseudocounts = pseudocounts

    @property
    def pwm(self):
        """Returns the positional weight matrix for the given collection."""
        return self._motif.pwm

    @property
    def sites(self):
        """Returns the binding sites in the collection."""
        return [str(instance) for instance in self._motif.instances]

    @property
    def site_count(self):
        """Returns the number of sites in the collection."""
        return len(self.sites)

    @property
    def length(self):
        """Returns the length of the sites."""
        return self._motif.length
