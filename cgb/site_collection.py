from Bio.Seq import Seq
from Bio import motifs
from Bio.motifs import jaspar
from Bio.Alphabet.IUPAC import unambiguous_dna


class SiteCollection:
    """The class definition for the collection of binding sites.

    SiteCollection encapsulates the transcription factor and its binding sites
    and constructs a position-independent probability model (PWM) which is used
    to build the TF binding model (see binding_model.py)

    """
    def __init__(self, sites, TF, pseudocounts=1):
        self._TF = TF
        instances = [Seq(site, unambiguous_dna) for site in sites]
        self._motif = motifs.create(instances)
        self._motif.pseudocounts = pseudocounts
        self._motif.name = self.TF.accession_number

    @property
    def TF(self):
        """Returns the TF object that the evidence is for."""
        return self._TF

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

    def to_jaspar(self, filename):
        """Writes the PWM to the given file in JASPAR format."""
        jaspar_motif = jaspar.Motif(matrix_id=self.TF.accession_number,
                                    name=self.TF.description,
                                    instances=self._motif.instances)
        with open(filename, 'w') as f:
            f.write(jaspar.write([jaspar_motif], 'jaspar'))
