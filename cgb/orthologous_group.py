"""Module for orthologous groups."""

import csv

from tqdm import tqdm

from my_logger import my_logger


class OrthologousGroup:
    """Class definition for OrthologousGroup.

    The OrthologousGroup class holds a group of genes which are orthologous to
    each other. Two genes are determined as orthologs if they are best BLAST
    hits for each other (reciprocal best BLAST hits).
    """
    def __init__(self, gene):
        self._genes = gene

    @property
    def genes(self):
        """Returns the list of orthologous genes."""
        return self._genes

    def member_from_genome(self, genome):
        """Returns the member of the group from the given genome.

        Returns None if the specified genome has no genes in the group.
        """
        genes = [g for g in self.genes if g.genome == genome]
        if genes:
            return genes[0]
        return None

    def __repr__(self):
        return str(self.genes)


# Class-associated functions
#
# The following functions provide the means to instantiate orthologous groups
# from a pre-defined subset of genes in all genomes under analysis and to
# export them in CSV format.


def construct_orthologous_groups(genes, genomes):
    """Constructs orthologous groups starting with the given list of genes.

    For each genome, candidate genes that are identified as likely to be
    regulated are tagged for orthology detection.

    This constructor function receives the genome objects and the list
    of genes from each of these genomes on which reciprocal BLAST will be
    applied to infer orhtologs.

    For each gene, it identifies the reciprocal best BLAST hits in other
    genomes and adds the gene and its orthologs to the orthologous group.
    Each orthologous group is a list of gene objects that have been found
    to be best-reciprocal BLAST hits.

    The function returns a list of orthologous groups
    """
    groups = []
    for gene in tqdm(genes):
        # Check whether gene is already in a group, if it is, it skips the gene
        # (continue goes back to for loop beginning)
        if any(gene in grp.genes for grp in groups):
            continue
        # If gene not in any group, create list of orthologous genes by
        # performing reciprocal BLAST against all genomes that are not the
        # gene's own genome
        rbhs = [gene.reciprocal_blast_hit(other_genome)
                for other_genome in genomes if gene.genome != other_genome]
        # Create the orthologous group with gene + orthologs on all other
        # genomes [if there are orthologs in the respective genomes]
        groups.append(OrthologousGroup([gene] + [rbh for rbh in rbhs if rbh]))
    return groups


def orthologous_groups_to_csv(groups, filename):
    genomes = list(set(g.genome for grp in groups for g in grp.genes))
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        header_row = [field for genome in genomes
                      for field in ['probability (%s)' % genome.strain_name,
                                    'locus_tag (%s)' % genome.strain_name,
                                    'product (%s)' % genome.strain_name]]
        csv_writer.writerow(header_row)
        for group in groups:
            row = []
            for genome in genomes:
                gene = group.member_from_genome(genome)
                if gene:
                    row.extend(['%.3f' % gene.operon.regulation_probability,
                                gene.locus_tag, gene.product])
                else:
                    row.extend(['', '', ''])
            csv_writer.writerow(row)
