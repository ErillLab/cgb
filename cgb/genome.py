import logging
import csv

from cached_property import cached_property
from tqdm import tqdm
from Bio.motifs import jaspar

from chromid import Chromid
from blast import BLAST
from pssm_model import PSSMModel
from misc import weighted_choice
from my_logger import my_logger


class Genome:
    def __init__(self, strain_name, accession_numbers):
        """Initializes the Genome object."""
        logging.debug('Creating genome: %s %s' %
                      (strain_name, str(accession_numbers)))
        self._strain_name = strain_name
        self._chromids = [Chromid(acc, self) for acc in accession_numbers]
        self._TF_instance = None
        self._TF_binding_model = None

    @property
    def strain_name(self):
        """Returns the name of the genome."""
        return self._strain_name

    @property
    def chromids(self):
        """Returns all chromosome/plasmid objects of the genome."""
        return self._chromids

    @property
    def num_chromids(self):
        """Returns the number of chromosome/plasmid objects."""
        return len(self.chromids)

    @cached_property
    def operons(self):
        """Returns all operons of the genome."""
        return [opr for chromid in self.chromids for opr in chromid.operons]

    def operons_to_csv(self, filename):
        """Writes all operons to the file in csv format."""
        with open(filename, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            header_row = ['chromid', 'start', 'end', 'strand', 'locus_tags',
                          'products']
            csv_writer.writerow(header_row)
            for opr in self.operons:
                row = [opr.chromid.accession_number,
                       opr.start,
                       opr.end,
                       opr.strand,
                       ', '.join(g.locus_tag for g in opr.genes),
                       ', '.join(g.locus_tag + ' (%s)' % g.product
                                 for g in opr.genes)]
                csv_writer.writerow(row)

    @cached_property
    def genes(self):
        """Returns all genes of the genome."""
        return [g for chromid in self.chromids for g in chromid.genes]

    @cached_property
    def protein_coding_genes(self):
        """Returns the protein coding genes of the genome."""
        return [g for chromid in self.chromids
                for g in chromid.protein_coding_genes]

    def genes_to_fasta(self):
        """Returns the sequences of all genes in FASTA format."""
        return '\n'.join(c.genes_to_fasta() for c in self.chromids)

    @cached_property
    def blast_client(self):
        """Returns the BLAST client.

        The target database is created with genes of the genome.
        """
        return BLAST(self.genes_to_fasta(), 'nucl', prefix=self.strain_name)

    @property
    def TF_instance(self):
        """Returns the instance of the TF in this genome.

        The TF-instance is identified through BLAST.
        """
        # TODO(sefa): make sure setter called before getter.
        return self._TF_instance

    @TF_instance.setter
    def TF_instance(self, value):
        """Sets the TF-instance attribute of the genome.

        Args:
            value (Protein): the transcription factor instance.
        """
        self._TF_instance = value

    @property
    def TF_binding_model(self):
        """Returns the TF-binding model inferred for this genome."""
        # TODO(sefa): make sure setter called before getter.
        return self._TF_binding_model

    def PSSM_model_to_jaspar(self, filename):
        """Writes the PSSM to the given file in jaspar format."""
        header = "genome: %s, TF: %s" % (self.strain_name, self.TF_instance.name)
        jaspar_motif = jaspar.Motif(matrix_id='', name=header,
                                    counts=self.TF_binding_model.pwm)
        with open(filename, 'w') as f:
            f.write(jaspar.write([jaspar_motif], 'jaspar'))

    def build_PSSM_model(self, collections, weights, prior_reg):
        """Builds a PSSM_model and sets the _TF_binding_model attribute.

        Args:
            collections ([SiteCollection]): list of site collections
            weights ([float]): list of weights, one per site collection
            prior_reg: prior probability of regulation of an operon
        Returns: None. Sets the _TF_binding_model attribute to the built model.
        """
        model = PSSMModel(collections, weights)
        random_seqs = self.random_seqs(length=model.length, count=100)
        bg_scores = [model.score_seq(random_seq) for random_seq in random_seqs]
        model.build_bayesian_estimator(bg_scores, prior_reg)
        self._TF_binding_model = model

    def random_seqs(self, length, count):
        """Returns random sequences drawn from the genome."""
        chromids = weighted_choice(self.chromids,
                                   weights=[c.length for c in self.chromids],
                                   count=count)
        return [c.random_seq(length=length) for c in chromids]

    def get_gene_by_locus_tag(self, locus_tag):
        """Returns the gene with the given locus tag."""
        gene, = [g for g in self.genes if g.locus_tag == locus_tag]
        return gene

    def find_gene_homolog(self, gene):
        """Returns the homolog gene of the genome."""
        blast_record = self.blast_client.tblastx(gene.to_fasta())
        locus_tag = self.blast_client.get_best_hit(blast_record)
        evalue = self.blast_client.get_e_value(blast_record)
        return self.get_gene_by_locus_tag(locus_tag), evalue

    def find_protein_homolog(self, protein):
        """Returns the homolog protein of the given protein."""
        # TODO(sefa): use blast database of protein coding genes only
        blast_record = self.blast_client.tblastn(protein.to_fasta())
        # TODO(sefa): exception when there are no hits
        locus_tag = self.blast_client.get_best_hit(blast_record)
        gene = self.get_gene_by_locus_tag(locus_tag)
        evalue = self.blast_client.get_e_value(blast_record)
        return gene.to_protein(), evalue

    def identify_TF_instance(self, proteins):
        """Finds the homolog of the given transcription factors.

        Args:
            proteins (list): List of proteins.
        Returns:
            protein: the homologous gene with the lowest evalue.
            None: if there are no homologs.
        """
        blast_hits = [self.find_protein_homolog(p) for p in proteins]
        if blast_hits:
            TF, _ = min(blast_hits, key=lambda x: x[1])
        else:
            my_logger.warning("No TF-instance found for %s. " %
                              self.strain_name)
        self._TF_instance = TF

    def scan_genome(self, threshold=0.5, filename=None):
        """Scans upstream regions of all operons for binding sites.

        Args:
            filename: csv filename to write results

        Returns:
            [(Operon, float)]: List of operons and their regulation
            probabilities, sorted by the probability.
        """
        scan_results = []
        for opr in tqdm(self.operons):
            p = opr.regulation_probability(self.TF_binding_model)
            if p >= threshold:
                scan_results.append((opr, p))
        scan_results.sort(key=lambda x: x[1], reverse=True)
        if filename:
            self._output_posterior_probabilities(scan_results, filename)

    def _output_posterior_probabilities(self, scan_results, filename):
        with open(filename, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            header_row = ['probability', 'operon_start', 'operon_end',
                          'operon_strand', 'genes', 'products']
            csv_writer.writerow(header_row)
            for opr, prob in scan_results:
                gene_names = ', '.join(g.locus_tag for g in opr.genes)
                products = ', '.join(g.locus_tag + ' (%s)' % g.product
                                     for g in opr.genes)
                row = ['%.3f' % prob, opr.start, opr.end, opr.strand, gene_names,
                       products]
                csv_writer.writerow(row)

    def __repr__(self):
        return (self.strain_name + ': ' +
                str([c.accession_number for c in self.chromids]))
