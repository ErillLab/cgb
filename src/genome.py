from chromid import Chromid


class Genome:
    def __init__(self, strain_name, accession_numbers):
        """Initializes the Genome object."""
        self._strain_name = strain_name
        self._chromids = [Chromid(acc, self) for acc in accession_numbers]

    @property
    def strain_name(self):
        """Returns the name of the genome."""
        return self._strain_name

    @property
    def chromids(self):
        """Returns all chromosome/plasmid objects of the genome."""
        return self._chromids

    @property
    def operons(self):
        """Returns all operons of the genome."""
        return [opr for chromid in self.chromids for opr in chromid.operons]

    @property
    def genes(self):
        """Returns all genes of the genome."""
        return [g for chromid in self.chromids for g in chromid.genes]

    def genes_to_fasta(self):
        """Returns the sequences of all genes in FASTA format."""
        return '\n'.join(c.genes_to_fasta() for c in self.chromids)

    def __repr__(self):
        return (self.strain_name + ': ' +
                str([c.accession_number for c in self.chromids]))
