import cStringIO
import random

from Bio import SeqIO
from cached_property import cached_property

import entrez_utils
import bio_utils
from gene import Gene
from operon import Operon
from misc import mean


class Chromid:
    """Definition for Chromid (chromosome//plasmid) class.

    Chromid class holds attributes and methods related to chromosome or plasmid
    of a genome.

    Given the accession number for the chromosome/plasmid, the corresponding
    NCBI GenBank record is fetched and parsed using the Biopython library. The
    record contains gene annotations (the location, locus tag, product, etc.)
    which are used to create Gene objects.

    The class implements a distance-based operon prediction method which is
    adaptive to the distribution of intergenic distances. See related method
    for more details.
    """

    def __init__(self, accession_number, genome):
        self._record = entrez_utils.get_genome_record(accession_number)
        self._genome = genome

    @property
    def genome(self):
        """Returns the Genome object that the chromid belongs to."""
        return self._genome

    @cached_property
    def record(self):
        """Returns the Biopython SeqRecord created from a GenBank file."""
        return SeqIO.read(cStringIO.StringIO(self._record), 'gb')

    @property
    def accession_number(self):
        """Returns the accession number."""
        return self.record.id

    @property
    def description(self):
        """Returns the description of the chromosome/plasmid."""
        return self.record.description

    @cached_property
    def sequence(self):
        """Returns the chromosome/plasmid sequence."""
        return str(self.record.seq)

    def random_seq(self, length):
        """Returns a random sequence drawn from the chromid sequence."""
        start = random.randint(0, self.length-length)
        return self.sequence[start:start+length]

    def random_seqs(self, length, count):
        """Returns random sequences drawn from the chromid."""
        return [self.random_seq(length) for _ in xrange(count)]

    @cached_property
    def length(self):
        """Returns the length of the genome sequence."""
        return len(self.sequence)

    def subsequence(self, start, end, strand=1):
        """Returns the specified DNA sequence."""
        seq = self.sequence[start:end]
        if strand == -1:
            seq = bio_utils.reverse_complement(seq)
        return seq

    @cached_property
    def genes(self):
        """Returns the list of genes of the chromosome/plasmid."""
        gene_list = []
        for f, next_f in zip(self.record.features, self.record.features[1:]):
            if f.type == 'gene':
                locus_tag = f.qualifiers['locus_tag']
                next_locus_tag = next_f.qualifiers.get('locus_tag')
                product_f = next_f if locus_tag == next_locus_tag else None
                gene_list.append(Gene(self, f, product_f))
        return gene_list

    @cached_property
    def protein_coding_genes(self):
        """Returns the protein coding genes of the chromosome/plasmid."""
        return [g for g in self.genes if g.product_type == 'CDS']

    @cached_property
    def operons(self):
        """Returns the list of operons of the chromosome/plasmid."""
        return self._operon_prediction()

    def genes_to_fasta(self):
        """Returns the sequences of all genes in FASTA format."""
        return '\n'.join(g.to_fasta() for g in self.genes)

    def _directons(self):
        """Returns the list of directons.

        A directon is a set of consecutive genes on one strand of DNA.
        """
        genes = sorted(self.genes, key=lambda g: g.start)
        directons = []
        cur_directon = [genes[0]]
        for cur_gene in genes[1:]:
            if cur_directon[-1].strand == cur_gene.strand:
                cur_directon.append(cur_gene)
            else:
                directons.append(cur_directon)
                cur_directon = [cur_gene]
        directons.append(cur_directon)
        return [directon if directon[0].is_forward_strand else directon[::-1]
                for directon in directons]

    def _operon_prediction(self):
        """Identifies all operons of the chromosome/plasmid.

        Two neighboring genes are in the same operon if their intergenic
        distance, which is defined as the mean intergenic distance between the
        first two genes of all opposite directons.
        """
        operons = []
        directons = self._directons()
        # Compute the mean intergenic distance of directons' first two genes.
        mean_dist = mean([directon[0].distance(directon[1])
                          for directon in directons if len(directon) > 1])

        directons_rest = self._directons()
        while directons_rest:
            processing = directons_rest
            directons_rest = []

            for directon in processing:
                operon = [directon[0]]
                i = 1
                while (i < len(directon) and
                       directon[i-1].distance(directon[i]) < mean_dist):
                    operon.append(directon[i])
                    i += 1
                operons.append(operon)
                if i < len(directon):
                    directons_rest.append(directon[i:])
        return [Operon(opr) for opr in operons]

    def find_closest_gene(self, pos):
        """Returns the closest gene and its distance to the pos.

        Negative distance means that the pos is downstream of the gene.
        """
        dist = (lambda gene: gene.start-pos if gene.strand == 1
                else pos-gene.end)
        dists = ((g, dist(g)) for g in self.genes)
        return min(dists, key=lambda x: abs(x[1]))

    def __repr__(self):
        return self.accession_number + ': ' + self.description
