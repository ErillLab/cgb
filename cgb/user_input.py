import json


class UserInput:
    def __init__(self, input_filename, config_filename=None):
        """Constructs Input object using the input file."""
        self._input = {}
        with open(input_filename) as f:
            for k, v in json.load(f).items():
                self._input[k] = v
        with open(config_filename) as f:
            self._input['config'] = {}
            for k, v in json.load(f).items():
                self._input['config'][k] = v

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
        directory = self._input['config']['log_dir']
        return directory
