import logging
import csv
from collections import namedtuple
import random

from cached_property import cached_property
from tqdm import tqdm
from Bio.motifs import jaspar

from .chromid import Chromid
from .blast import BLAST, BlastNoHitFoundException
from .pssm_model import PSSMModel
from .misc import weighted_choice
from .misc import temp_file_name
from .misc import mean
from .my_logger import my_logger
from .bio_utils import reverse_complement
from .bio_utils import weblogo
from .gene import NotProteinCodingGeneException


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

        logging.info('Creating genome: %s %s' %
                     (strain_name, str(accession_numbers)))
        self._strain_name = strain_name
        self._chromids = [Chromid(acc, self) for acc in accession_numbers]
        self._TF_instance = None   # TF-instance in this genome
        self._TF_binding_model = None  # binding model tailored for the genome

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
    def directons(self):
        """Returns the list of directons in all chromids."""
        return [d for chromid in self.chromids for d in chromid.directons]

    def intergenic_distance_threshold(self, sigma=1.0):
        """The intergenic distance threshold used for operon prediction.

        Two adjacent genes on the same strand belong to the same operon if
        the distance between them is less than the determined threshold.

        Args: sigma (float): scale factor to calibrate the operon prediction
            distance threshold.
        """
        dists = [d[0].distance(d[1]) for d in self.directons if len(d) > 1]
        return sigma * mean(dists)

    @property
    def operons(self):
        """Returns all operons of the genome."""
        oprs = [opr for chromid in self.chromids for opr in chromid.operons]
        for i, opr in enumerate(oprs, start=1):
            opr.set_operon_id = i
        return oprs

    def remove_operons(self):
        """Removes predicted operons for each chromid.
        """
        for chromid in self.chromids:
            #if at least one operon has been predicted in chromid
            if len(chromid.operons)>0:
                chromid.remove_operons()
                del self.operons[:]

    def operon_prediction(self, *args):
        """Predicts operons for each chromid.

        See chromid.operon_prediction for details.
        """
        start_id = 1            # The id for the first operon
        for chromid in self.chromids:
            chromid.operon_prediction(start_id, *args)
            #if at least one operon has been predicted in chromid
            if len(chromid.operons)>0:
                start_id = chromid.operons[-1].operon_id + 1

            my_logger.info("Number of operons (%s): %d" %
                          (self.strain_name, self.num_operons))


    @cached_property
    def length(self):
        """Returns the total length of the genome."""
        return sum(c.length for c in self.chromids)

    @property
    def num_operons(self):
        """Returns the number of operons of the genome."""
        return len(self.operons)

    def operons_to_csv(self, filename):
        """Writes all operons to the file in csv format.

        Args:
            filename (string): the CSV filename to write all operons to
        """
        with open(filename, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            header_row = ['operon_id', 'chromid', 'start', 'end', 'strand',
                          'locus_tags', 'products']
            csv_writer.writerow(header_row)
            for opr in self.operons:
                row = [opr.operon_id,
                       opr.chromid.accession_number,
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

        Uses the list of genes in the current genome to invoke the BLAST
        object constructor, which will call BLAST to create a genome-specific
        database and provides access to BLAST methods.
        """
        return BLAST(self.genes_to_fasta(), 'nucl', prefix=self.strain_name)

    @property
    def TF_instance(self):
        """Returns the instance of the TF in this genome.

        The TF-instance should be identified through BLAST before calling.
        """
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
        """Returns the TF-binding model inferred for this genome.

        The TF-binding model tailored for this genome is set via
        "build_PSSM_model" method.
        """
        return self._TF_binding_model

    def build_PSSM_model(self, collections, weights):
        """Builds a PSSM_model and sets the _TF_binding_model attribute.

        Args:
            collections ([SiteCollection]): list of site collections
            weights ([float]): list of weights, one per site collection
        Returns: None. Sets the _TF_binding_model attribute to the built model.
        """
        my_logger.info("Building PSSM model (%s)" % self.strain_name)
        # Create a PSSM model using site collections and associated weights.
        model = PSSMModel(collections, weights)
        # Randomly select sites from promoter regions
        noncoding_regions = [g.upstream_noncoding_region_sequence()
                             for g in self.genes]
        seqs = [seq[start:start+model.length] for seq in noncoding_regions
                for start in xrange(len(seq)-model.length)]
        sampled_seqs = random.sample(seqs, min(10000, len(seqs)))
        # Score randomly selected sites to estimate the bg distribution
        bg_scores = []
        for seq in sampled_seqs:
            bg_scores.append(model.score_seq(seq))
        # Create a Bayesian estimator which is used to compute the probability
        # of TF-binding of any given sequence.
        model.build_bayesian_estimator(bg_scores)
        self._TF_binding_model = model

    def random_seqs(self, length, count):
        """Returns random sequences drawn from the genome."""
        # Choose the chromosomes randomly, weighted by the length.
        chromids = weighted_choice(self.chromids,
                                   weights=[c.length for c in self.chromids],
                                   count=count)
        return [c.random_seq(length=length) for c in chromids]

    def get_gene_by_locus_tag(self, locus_tag):
        """Returns the gene with the given locus tag.

        Args:
            locus_tag (string): the locus_tag of the gene.
        """
        gene, = [g for g in self.genes if g.locus_tag == locus_tag]
        return gene

    def find_gene_homolog(self, gene):
        """Returns the Gene object that is homologous to the given gene.

        Invokes TBLASTX to identify the best hit of the query gene in the
        genome. Requires the BLAST package to be installed and that the BLAST
        package binaries are in the path.

        Args:
            gene (Gene): the query gene.
        Returns:
            (Gene, float): The best BLAST hit in the genome
            for the query and its BLAST e-value.
        """
        # Perform a tblastx search with the given gene against the genome.
        # The blast_client returns a Biopython blast_record.
        blast_record = self.blast_client.tblastx(gene.to_fasta())
        # Call the get_best_hit method to get the locus_tag and e-value of the
        # first BLAST hit.
        locus_tag = self.blast_client.get_best_hit(blast_record)
        evalue = self.blast_client.get_e_value(blast_record)
        return self.get_gene_by_locus_tag(locus_tag), evalue

    def find_protein_homolog(self, protein):
        """Returns the homolog protein of the query protein in this genome.

        Invokes TBLASTN to identify the best hit of the query protein in the
        genome. Requires the BLAST package to be installed and that the BLAST
        binaries are in the path.

        Args:
            protein (Protein): the query protein.
        Returns:
            (Protein, float): The homologous protein and the BLAST e-value.
        """
        # Perform a tblastn search with the given protein against the genome.
        # The blast_client returns a Biopython blast_record.
        blast_record = self.blast_client.tblastn(protein.to_fasta())
        # Get the best hit to get the locus_tag and e-value of the first
        # BLAST hit.
        locus_tag = self.blast_client.get_best_hit(blast_record)
        gene = self.get_gene_by_locus_tag(locus_tag)
        evalue = self.blast_client.get_e_value(blast_record)
        protein = gene.to_protein()
        return protein, evalue

    def identify_TF_instance(self, proteins):
        """Finds the instance of the transcription factor of interest.

        Given a list of proteins corresponding to the TF of interest (provided
        by the user; the list of proteins for which binding motifs are
        provided), it calls "find_protein_homolog" method to identify the
        homologous protein in this genome. Among the returned BLAST hits, it
        picks the lowest e-value result to define the orthologous TF in the
        genome. If the lowest e-value does not meet the threshold, this is
        reported to the user and the genome is excluded from analysis.

        Args:
            proteins (list): List of TF-instances provided by the user.
        Returns:
            protein: the homologous gene with the lowest evalue.
            None: if there are no homologs.

        """
        # Find best BLAST hit for each instance of the TF provided by user (p)
        blast_hits = []
        for p in proteins:
            # Identify the protein of the genome that is homologous to p.
            try:
                homolog_protein = self.find_protein_homolog(p)
                blast_hits.append(homolog_protein)
            except (BlastNoHitFoundException, NotProteinCodingGeneException):
                pass
        if blast_hits:
            # If there are BLAST hits, return the one with the best e-value.
            TF, _ = min(blast_hits, key=lambda x: x[1])
            my_logger.info("%s" % TF.accession_number)
        else:
            # Otherwise, set the TF-instance to None.
            # The genome will be dropped from the analysis.
            TF = None
            my_logger.warning("No TF-instance found for %s. " %
                              self.strain_name)

        self._TF_instance = TF

    def calculate_regulation_probabilities(self, prior,user_input):
        """Calculates posterior probability of regulation for _all_ genes."""
        for gene in tqdm(self.genes):
            gene.calculate_regulation_probability(prior,user_input)

    def infer_regulons(self, threshold=0.5, filename=None):
        """Scans upstream regions of all operons for binding sites.

        Args:
            threshold (float): The threshold for posterior probability of
                binding. Only the operons with a TF-binding probability higher
                than the threshold are reported.
            filename (sting): The CSV file to write the computed posterior
                probabilities, one line per operon.

        Returns:
            [(Operon, float)]: List of operons and their regulation
            probabilities, sorted by the probability.
        """
        my_logger.info("Identifying putative regulons (%s)" % self.strain_name)
        # Find regulated operons
        results = [(opr, opr.regulation_probability) for opr in self.operons]
        # Sort results by posterior probability.
        results.sort(key=lambda x: x[1], reverse=True)
        # Output results to the CSV file, if filename is provided.
        if filename:
            self._output_posterior_probabilities(results, filename)
        # filter regulons by the probability threshold
        return [opr for (opr, p) in results if p >= threshold]

    def _output_posterior_probabilities(self, results, filename):
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
            for opr, prob in results:
                gene_names = ', '.join(g.locus_tag for g in opr.genes)
                products = ', '.join(g.locus_tag + ' (%s)' % g.product
                                     for g in opr.genes)
                row = ['%.3f' % prob, opr.start, opr.end, opr.strand,
                       gene_names, products]
                csv_writer.writerow(row)

    @property
    def putative_sites(self):
        """Returns the list of putative sites in non-coding regions."""
        # 'identify_sites' method should be called before.
        return self._putative_sites

    @property
    def weblogo_from_putative_sites(self):
        """Returns the file of sequence logo built using putative sites."""
        filename = temp_file_name()
        weblogo([site.chromid.subsequence(site.start, site.end, site.strand)
                 for site in self.putative_sites], filename)
        return filename

    def identify_sites(self, user_input, filename=None):
        """Returns the list of sites in non-coding regions.

        It searches exclusively the [-promoter_up, promoter_dw] for all genes
        in the genome. It returns all sites with a score over threshold, tags
        found sites as operator or intergenic (if distance is >300, or whatever
        the user specifies). Finally, it reports identified sites to a CSV file
        if filename is provided.

        The reason for identifying binding sites upstream of all genes, rather
        than only for operons, is that the binding site predictions are used to
        improve the operon prediction. That is, if there is a putative binding
        site upstream of an in-between gene in an operon, the operon is split
        into two by the gene.

        Args:
            user_input: to get the promoter_up and promoter_dw max-distances
            from the start of the gene to report binding sites.

            filename (string): the CSV file to report putative binding sites.
        """
        threshold = self.TF_binding_model.threshold()  # score threshold
        site_len = self.TF_binding_model.length  # Length of the binding sites
        sites = []
        my_logger.info("Identifying sites in %s" % self.strain_name)
        for gene in tqdm(self.genes):
            # Locate the upstream non-coding region.
            # Provide None as promoter_up_distance, so that the whole intergenic
            # region up the next gene TLS is scanned
            start, end = gene.upstream_noncoding_region_location(None,
                                                user_input.promoter_dw_distance)
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

        # Sort the identified sites by their scores.
        sites.sort(key=lambda site: site.score, reverse=True)
        # Store the results in a private attribute.
        self._putative_sites = sites

        if filename:
            self._output_identified_sites(sites, 
                                      user_input.promoter_up_distance, filename)

    def _output_identified_sites(self, sites, promoter_up, filename):
        """Reports the idenitied sites to a CSV file.
        Args:
            sites ([Site]): List of Site tuples to be reported
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
        header = self.strain_name
        jaspar_motif = jaspar.Motif(matrix_id='', name=header,
                                    counts=self.TF_binding_model.pwm)
        with open(filename, 'w') as f:
            f.write(jaspar.write([jaspar_motif], 'jaspar'))

    def __repr__(self):
        return (self.strain_name + ': ' +
                str([c.accession_number for c in self.chromids]))
