import logging
import csv
from collections import namedtuple

from cached_property import cached_property
from tqdm import tqdm
from Bio.motifs import jaspar

from chromid import Chromid
from blast import BLAST
from pssm_model import PSSMModel
from misc import weighted_choice
from my_logger import my_logger
from bio_utils import reverse_complement

Site = namedtuple('Site', 'chromid start end strand score gene')


class Genome:
    """Definition for Genome class.

    Genome class encapsulates all chromosomes and plasmids of a species.

    In addition to providing access to all chromid attributes and methods, it
    holds a BLAST client to query the database constructed with all the genes
    in the genome.

    The class also contains a genome-specific TF-binding model, which is used
    to identify individual binding sites on the genome, as well as to infer the
    posterior probabilities of regulation for each operon. It also provides
    methods to output predicted binding sites and posterior probability
    estimations to corresponding comma separated value (CSV) files.
    """
    def __init__(self, strain_name, accession_numbers):
        """Initializes the Genome object.

        Args:
            strain_name (string): User-provided name for the strain
            accession_numbers (list): List of corresponding accession numbers.
        """

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
        """Returns the homolog gene of the genome.

        Args:
            gene (Gene): the query gene.
        Returns:
            (Gene, float): The homologous gene and the BLAST e-value.
        """
        blast_record = self.blast_client.tblastx(gene.to_fasta())
        locus_tag = self.blast_client.get_best_hit(blast_record)
        evalue = self.blast_client.get_e_value(blast_record)
        return self.get_gene_by_locus_tag(locus_tag), evalue

    def find_protein_homolog(self, protein):
        """Returns the homolog protein of the given protein.

        Args:
            protein (Protein): the query protein.
        Returns:
            (Protein, float): The homologous protein and the BLAST e-value.
        """
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

    def infer_regulation(self, threshold=0.5, filename=None):
        """Scans upstream regions of all operons for binding sites.

        Args:
            filename: csv filename to write posterior probabilities

        Returns:
            [(Operon, float)]: List of operons and their regulation
            probabilities, sorted by the probability.
        """
        # Find regulated operons
        regulons = []
        for opr in tqdm(self.operons):
            p = opr.regulation_probability(self.TF_binding_model)
            if p >= threshold:
                regulons.append((opr, p))
        regulons.sort(key=lambda x: x[1], reverse=True)
        if filename:
            self._output_posterior_probabilities(regulons, filename)
        return regulons

    def _output_posterior_probabilities(self, scan_results, filename):
        """Outputs posterior probabilities to a csv file.

        Args:
            scan_results (list): List of (Operon, probability) pairs.
            filename (string): the path to the output CSV file.
        """
        with open(filename, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            header_row = ['probability', 'operon_start', 'operon_end',
                          'operon_strand', 'genes', 'products']
            csv_writer.writerow(header_row)
            for opr, prob in scan_results:
                gene_names = ', '.join(g.locus_tag for g in opr.genes)
                products = ', '.join(g.locus_tag + ' (%s)' % g.product
                                     for g in opr.genes)
                row = ['%.3f' % prob, opr.start, opr.end, opr.strand,
                       gene_names, products]
                csv_writer.writerow(row)

    @property
    def putative_sites(self):
        """Returns the lsit of putative sites in non-coding regions."""
        return self._putative_sites

    def identify_sites(self, promoter_up=300, filename=None):
        """Returns the list of sites in non-coding regions..

        It searches exclusively the [up to previous gene, +50] for all genes in
        the genome. It returns all sites with a score over threshold, and tags
        the found sites as operator or intergenic (if distance is >300, or
        whatever the user specifies). Finally, it reports identified sites to a
        CSV file if filename is provided.
        """
        if hasattr(self, '_putative_sites' ):
            return self._putative_sites

        threshold = self.TF_binding_model.threshold()  # score threshold
        site_len = self.TF_binding_model.length  # Length of the binding sites
        sites = []
        my_logger.debug("Identifying sites in %s" % self.strain_name)
        for gene in tqdm(self.genes):
            start, end = gene.upstream_noncoding_region()
            # Score forward strand
            seq = gene.chromid.subsequence(start, end)
            scores = self.TF_binding_model.score_seq(seq, both=False)
            for i, score in enumerate(scores, start=start):
                if score >= threshold:
                    sites.append(
                        Site(gene.chromid, i, i+site_len, 1, score, gene))
            # Score reverse strand
            rc_seq = reverse_complement(seq)
            rc_scores = self.TF_binding_model.score_seq(rc_seq, both=False)
            for i, score in enumerate(reversed(rc_scores), start=start):
                if score >= threshold:
                    sites.append(
                        Site(gene.chromid, i, i+site_len, -1, score, gene))
        sites.sort(key=lambda site: site.score, reverse=True)
        self._putative_sites = sites

        if filename:
            self._output_identified_sites(sites, promoter_up, filename)

    def _output_identified_sites(self, sites, promoter_up, filename):
        """Reports the idenitied sites to a CSV file.
        Args:
            sites (list(Site)): List of Site tuples to be reported
            promoter_up (float): the distance threshold for operator region. If
            the distance to the gene is <threshold, the site is in operator
            region. Otherwise, it is in intergenic region.
            filename (string): the path to csv file.
        """
        with open(filename, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            header_row = ['score', 'site', 'chromid', 'start', 'end', 'strand',
                          'distance', 'category', 'gene_strand', 'gene_start',
                          'gene_end', 'gene_locus_tag', 'gene_product']
            csv_writer.writerow(header_row)
            for site in sites:
                dist = site.gene.distance_to_region(site.start, site.end)
                category = 'operator' if dist < promoter_up else 'intergenic'
                site_seq = site.chromid.subsequence(
                    site.start, site.end, site.strand)
                row = [site.score, site_seq, site.chromid.accession_number,
                       site.start, site.end, site.strand, dist, category,
                       site.gene.strand, site.gene.start,
                       site.gene.end, site.gene.locus_tag,
                       site.gene.product]
                csv_writer.writerow(row)

    def output_TF_binding_model(self, filename):
        """Writes the model to the given file."""
        self._PSSM_model_to_jaspar(filename)

    def _PSSM_model_to_jaspar(self, filename):
        """Writes the PSSM to the given file in jaspar format."""
        header = "genome: %s, TF: %s" % (self.strain_name, self.TF_instance.name)
        jaspar_motif = jaspar.Motif(matrix_id='', name=header,
                                    counts=self.TF_binding_model.pwm)
        with open(filename, 'w') as f:
            f.write(jaspar.write([jaspar_motif], 'jaspar'))

    def __repr__(self):
        return (self.strain_name + ': ' +
                str([c.accession_number for c in self.chromids]))
