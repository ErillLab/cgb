from cached_property import cached_property

from chromid import Chromid
from blast import BLAST


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

    @property
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
        return BLAST(self.genes_to_fasta(), 'nucl')

    def get_gene_by_locus_tag(self, locus_tag):
        """Returns the gene with the given locus tag."""
        gene, = [g for g in self.genes if g.locus_tag == locus_tag]
        return gene

    def find_gene_homolog(self, gene):
        """Returns the homolog gene of the genome."""
        blast_record = self.blast_client.tblastx(gene.to_fasta())
        # TODO(sefa): exception when there are no hits
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
        return min(blast_hits, key=lambda x: x[1]) if blast_hits else None

    def __repr__(self):
        return (self.strain_name + ': ' +
                str([c.accession_number for c in self.chromids]))
