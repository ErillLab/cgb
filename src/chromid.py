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
    """Class for chromosomes and plasmids."""
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

    def random_seqs(self, length, count):
        """Returns random sequences drawn from the genome."""
        starts = [random.randint(0, self.length-length) for _ in xrange(count)]
        return [self.sequence[start:start+length] for start in starts]

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
        return [Gene(f, self) for f in self.record.features if f.type=='gene']

    @cached_property
    def operons(self):
        """Returns the list of operons of the chromosome/plasmid."""
        return self._operon_prediction()

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

    def __repr__(self):
        return self.accession_number + ': ' + self.description
