import os

from genome import Genome
from protein import Protein
from site_collection import SiteCollection
from my_logger import my_logger
from phylo import Phylo
from user_input import UserInput
from orthologous_group import construct_orthologous_groups
from orthologous_group import orthologous_groups_to_csv


def directory(*paths):
    d = os.path.join(*paths)
    if not os.path.exists(d):
        os.makedirs(d)
    return d


def create_genomes(user_input):
    my_logger.info("Started: create genomes")
    genomes = [Genome(name, accessions)
               for name, accessions in user_input.genome_name_and_accessions]
    # write operons to files
    log_dir = directory(user_input.log_dir, 'operons')
    for g in genomes:
        g.operons_to_csv(os.path.join(log_dir, g.strain_name+'_operons.csv'))
    my_logger.info("Finished: create genomes")
    return genomes


def identify_TF_instance_in_genomes(genomes, proteins):
    """Identifies TF instance for each genome based on given proteins."""
    for g in genomes:
        g.identify_TF_instance(proteins)


def set_TF_binding_models(user_input, genomes, site_collections, weights):
    """Sets the TF-binding model for each genome."""
    prior_reg = 0.01            # TODO(sefa): get it from input file
    for g in genomes:
        g.build_PSSM_model(site_collections, weights, prior_reg)
        log_dir = directory(user_input.log_dir, 'derived_PSWM')
        g.output_TF_binding_model(os.path.join(log_dir, g.strain_name+'.jaspar'))


def create_proteins(user_input):
    my_logger.info("Started: create proteins")
    proteins = [Protein(accession, user_input.TF_name)
                for accession in user_input.protein_accessions]
    my_logger.info("Finished: create proteins")
    return proteins


def create_site_collections(user_input, proteins):
    my_logger.info("Started: create site collections")
    collections = [SiteCollection(sites, protein)
                   for sites, protein in zip(user_input.sites_list, proteins)]
    my_logger.debug("Writing site collections.")
    log_dir = directory(user_input.log_dir, 'user_PSWM')
    for collection in collections:
        collection.to_jaspar(
            os.path.join(log_dir, collection.TF.accession_number+'.jaspar'))
    my_logger.info("Finished: create site collections")
    return collections


def simple_motif_weighting(site_collections):
    """Weights collections according to collection sizes."""
    return [collection.site_count for collection in site_collections]


def phylogenetic_weighting(site_collections, phylogeny):
    # TODO(sefa): implement this.
    return uniform_weighting(site_collections)


def uniform_weighting(site_collections):
    return [1 for _ in site_collections]


def compute_motif_weights(site_collections, weighting_scheme):
    my_logger.debug("TF-binding evidence weighting scheme: %s" %
                    weighting_scheme)
    assert weighting_scheme in ['simple', 'phylogenetic']
    if weighting_scheme == 'simple':
        return simple_motif_weighting(site_collections)
    elif weighting_scheme == 'phylogenetic':
        return phylogenetic_weighting(site_collections)
    else:
        return uniform_weighting(site_collections)


def genome_scan_results_to_genes(genome_scan):
    """Returns the list of regulated genes."""
    return [g for opr, _ in genome_scan for g in opr.genes]


def infer_regulations(user_input, genomes):
    """Scans genomes for TF-binding sites and writes results to csv files."""
    log_dir = directory(user_input.log_dir, 'posterior_probs')
    threshold = 0.5             # TODO(sefa): read this from input file
    all_regulated_genes = []
    for g in genomes:
        report_filename = os.path.join(log_dir, g.strain_name+'.csv')
        genome_scan = g.infer_regulation(threshold, report_filename)
        all_regulated_genes.extend(genome_scan_results_to_genes(genome_scan))
    return all_regulated_genes


def search_sites(user_input, genomes):
    log_dir = directory(user_input.log_dir, 'identified_sites')
    for g in genomes:
        my_logger.debug("threshold: %.2f" % g.TF_binding_model.threshold())
        report_filename = os.path.join(log_dir, g.strain_name+'.csv')
        sites = g.search_sites(report_filename)
    return sites


def create_orthologous_groups(user_input, genes, genomes):
    groups = construct_orthologous_groups(genes, genomes)
    # Write groups to file
    log_dir = directory(user_input.log_dir)
    orthologous_groups_to_csv(groups, os.path.join(log_dir, 'orthologs.csv'))
    return groups


def create_phylogeny(genomes, proteins, user_input):
    phylo = Phylo(proteins + [g.TF_instance for g in genomes])
    phylo.to_newick(os.path.join(user_input.log_dir, 'phylogeny.nwk'))
    return phylo


def main():
    # Read user input
    user_input = UserInput('../tests/input.json', '../tests/config.json')
    # Create proteins
    proteins = create_proteins(user_input)
    # Create genomes
    genomes = create_genomes(user_input)
    # Identify TF instances for each genome.
    identify_TF_instance_in_genomes(genomes, proteins)
    # Create binding evidence
    site_collections = create_site_collections(user_input, proteins)
    collection_weights = compute_motif_weights(site_collections, 'simple')
    # Set TF-binding model for each genome.
    set_TF_binding_models(user_input, genomes, site_collections,
                          collection_weights)
    # Score genomes
    all_regulated_genes = infer_regulations(user_input, genomes)
    sites = search_sites(user_input, genomes)
    # Create orthologous groups
    grps = create_orthologous_groups(user_input, all_regulated_genes, genomes)
    # Create phylogeny
    phylo = create_phylogeny(genomes, proteins, user_input)
