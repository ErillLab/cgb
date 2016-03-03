import json


class UserInput:
    def __init__(self, filename):
        """Constructs Input object using the input file.

        Args:
            filename (string): path to a file in JSON format.
        """
        with open(filename) as f:
            self._input = json.load(f)

    @property
    def genome_name_and_accessions(self):
        """Returns the list of (genome name, accession_numbers) tuples."""
        return [(g['name'], g['accession_numbers'])
                for g in self._input['genomes']]

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
    def log_dir(self):
        """Returns the directory to be used for logging."""
        directory = self._input['configuration']['log_dir']
        return directory
