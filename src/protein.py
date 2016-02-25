import cStringIO

from cached_property import cached_property
from Bio import SeqIO

import entrez_utils


class Protein:
    def __init__(self, accession_number):
        self._record = entrez_utils.get_protein_record(accession_number)

    @cached_property
    def record(self):
        """Returns the Biopython SeqRecord created from the NCBI record."""
        return SeqIO.read(cStringIO.StringIO(self._record), 'gb')

    @property
    def accession_number(self):
        """Returns the accession number of the protein."""
        return self.record.id

    @property
    def description(self):
        """Returns the description of the protein."""
        return self.record.description

    @property
    def sequence(self):
        """Returns the amino acid sequence of the protein."""
        return str(self.record.seq)
