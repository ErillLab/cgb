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

    def add_to_group(self, gene):
        self._genes.append(gene)

    def membership_test(self, other):
        return any(other.reciprocal_blast_hit(g.genome) == g
                   for g in self.genes)

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


def construct_orthologous_groups(genes, genomes):
    """Constructs orthologous groups starting with the given list of genes.

    For each gene, it identifies the reciprocal best BLAST hits in other
    genomes and adds them to the orthologous groups as well.
    """
    groups = []
    for gene in tqdm(genes):
        if any(gene in grp.genes for grp in groups):
            continue
        rbhs = [gene.reciprocal_blast_hit(other_genome)
                for other_genome in genomes if gene.genome != other_genome]
        groups.append(OrthologousGroup([gene] + [rbh for rbh in rbhs if rbh]))
    return groups


def orthologous_groups_to_csv(groups, filename):
    genomes = list(set(g.genome for grp in groups for g in grp.genes))
    with open(filename, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        header_row = [field for genome in genomes
                      for field in ['locus_tag (%s)' % genome.strain_name,
                                    'product (%s)' % genome.strain_name]]
        csv_writer.writerow(header_row)
        for group in groups:
            row = []
            for genome in genomes:
                gene = group.member_from_genome(genome)
                if gene:
                    row.extend([gene.locus_tag, gene.product])
                else:
                    row.extend(['', ''])
            csv_writer.writerow(row)
