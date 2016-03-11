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
    """Makes a directory.

    Creates the directory if it doesn't exist.

    Args:
        paths (list(string)): the full path for the directory
    Returns:
        string: the full path for the directory
    """
    d = os.path.join(*paths)
    if not os.path.exists(d):
        os.makedirs(d)
    return d


def create_genomes(user_input):
    """Creates all genomes used in the analysis.

    For each genome, it writes operons to csv files.
    """
    my_logger.info("Started: create genomes")
    genomes = [Genome(name, accessions)
               for name, accessions in user_input.genome_name_and_accessions]
    my_logger.info("Finished: create genomes")
    return genomes


def output_operons(user_input, genome):
    # write operons to files
    log_dir = directory(user_input.log_dir, 'operons')
    genome.operons_to_csv(
        os.path.join(log_dir, genome.strain_name+'_operons.csv'))


def identify_TF_instance_in_genomes(genomes, proteins):
    """Identifies TF instance for each genome based on given proteins."""
    for g in genomes:
        my_logger.info("Idenifying TF instance for /%s/" % g.strain_name)
        g.identify_TF_instance(proteins)


def set_TF_binding_model(user_input, genome, site_collections, weights):
    """Sets the TF-binding model for each genome.

    The TF_binding_models are genome-specific and built using the collections
    of sites and their associated weights.

    For each genome, built models are written to corresponding CSV files.

    Args:
        user_input (UserInput): options provided by the user
        genome (Genome): genome of interest
        site_collections (list(SiteCollection)): list of collections
        weights (list(float)): weights used for combining evidence
    """
    genome.build_PSSM_model(site_collections, weights)
    log_dir = directory(user_input.log_dir, 'derived_PSWM')
    genome.output_TF_binding_model(
        os.path.join(log_dir, genome.strain_name+'.jaspar'))


def create_proteins(user_input):
    """Creates all proteins used in the analysis."""
    my_logger.info("Started: create proteins")
    proteins = [Protein(accession, user_input.TF_name)
                for accession in user_input.protein_accessions]
    my_logger.info("Finished: create proteins")
    return proteins


def create_site_collections(user_input, proteins):
    """Creates SiteCollection objects from the user-provided collections."""
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


def phylogenetic_weighting(site_collections, genome, phylogeny):
    dists = [phylogeny.distance(collection.TF, genome.TF_instance)
             for collection in site_collections]
    normalized_dists = [float(d)/sum(dists) for d in dists]
    similarities = [1-nd for nd in normalized_dists]
    weights = [s/sum(similarities) for s in similarities]
    return weights


def uniform_weighting(site_collections):
    return [1 for _ in site_collections]


def compute_motif_weights(genome, site_collections, weighting_scheme,
                          phylogeny):
    my_logger.info("TF-binding evidence weighting scheme: %s" %
                   weighting_scheme)
    assert weighting_scheme in ['simple', 'phylogenetic']
    # Compute genome-specific weights of each site collection.
    if weighting_scheme == 'simple':
        weights = simple_motif_weighting(site_collections)
    elif weighting_scheme == 'phylogenetic':
        weights = phylogenetic_weighting(site_collections, genome, phylogeny)
    else:
        weights = uniform_weighting(site_collections)
    # Normalize weights
    weights = [float(w)/sum(weights) for w in weights]
    my_logger.info("Binding evidence weights for /%s/: %s" %
                   (genome.strain_name, str(weights)))
    return weights


def genome_scan_results_to_genes(genome_scan):
    """Returns the list of regulated genes."""
    return [g for opr, _ in genome_scan for g in opr.genes]


def get_prior(genome, user_input, weights):
    """Infers the prior probability of binding.

    This is done by predicting operons in the genome of the TF for which the
    motif is provided and dividing the number of sites in the motif (
    assuming that there is one site per operon) by the number of operons. For
    instance, if 30 sites are available for LexA in E. coli, then the prior for
    regulation is 30/~2300 operons.

    Then it is assigned to genome using the same strategy defined for mixing
    the TF-binding models. That is, simple combination (arithmetic mean of
    priors) or weighted combination (with phylogenetic distances).
    """
    if user_input.prior_regulation_probability:
        return user_input.prior_regulation_probability

    my_logger.info("Prior probability not provided, "
                   "inferring from the provided motifs. /%s"
                   % genome.strain_name)
    site_collections = genome.TF_binding_model.site_collections
    priors = [float(collection.site_count)/genome.num_operons
              for collection in site_collections]
    return sum(p*w for p, w in zip(priors, weights))


def infer_regulations(user_input, genome, prior):
    """Scans genomes for TF-binding sites and writes results to csv files."""
    log_dir = directory(user_input.log_dir, 'posterior_probs')
    threshold = user_input.probability_threshold
    my_logger.info("Regulation probability threshold: %.2f" % threshold)
    my_logger.info("Inferring regulations for /%s/." % genome.strain_name)
    my_logger.info("Prior probability of regulation: %f" % prior)
    report_filename = os.path.join(log_dir, genome.strain_name+'.csv')
    genome_scan = genome.infer_regulation(prior, threshold, report_filename)
    return genome_scan_results_to_genes(genome_scan)


def search_sites(user_input, genome):
    """Searches for putative binding sites.

    It uses genome's TF-binding model for binding site search and outputs
    the results into a CSV file accordingly.
    """
    log_dir = directory(user_input.log_dir, 'identified_sites')
    my_logger.info("Scoring genome /%s/." % genome.strain_name)
    my_logger.info(
        "Site score threshold: %.2f" % genome.TF_binding_model.threshold())
    report_filename = os.path.join(log_dir, genome.strain_name+'.csv')
    sites = genome.identify_sites(filename=report_filename)
    return sites


def create_orthologous_groups(user_input, genes, genomes):
    """Creates orthologous groups from the given list of genes.

    For each gene in the list, it searches for orthologs in other genomes and
    outputs results to the CSV file.
    """
    my_logger.info("Creating orthologous groups")
    groups = construct_orthologous_groups(genes, genomes)
    # Write groups to file
    log_dir = directory(user_input.log_dir)
    orthologous_groups_to_csv(groups, os.path.join(log_dir, 'orthologs.csv'))
    return groups


def create_phylogeny(genomes, proteins, user_input):
    """Creates a phylogeny from the given proteins.

    In addition to transcription factors with binding evidence, their
    homologous counterparts in the analyzed genomes are used to build the
    phylogeny.
    """
    phylo = Phylo(proteins + [g.TF_instance for g in genomes])
    phylo.to_newick(os.path.join(user_input.log_dir, 'phylogeny.nwk'))
    phylo.draw_ascii()
    return phylo


def main():
    """The entry-point for the pipeline."""
    # Read user input
    user_input = UserInput('../tests/input.json', '../tests/config.json')

    # Create proteins
    proteins = create_proteins(user_input)

    # Create genomes
    genomes = create_genomes(user_input)

    # Identify TF instances for each genome.
    identify_TF_instance_in_genomes(genomes, proteins)

    # Create phylogeny
    phylogeny = create_phylogeny(genomes, proteins, user_input)

    # Create binding evidence
    site_collections = create_site_collections(user_input, proteins)

    all_regulated_genes = []
    for genome in genomes:
        # Create weights for combining motifs and inferring priors, if not set
        weights = compute_motif_weights(genome, site_collections,
                                        user_input.weighting_scheme, phylogeny)
        # Set TF-binding model for each genome.
        set_TF_binding_model(user_input, genome, site_collections, weights)
        # Score genome
        search_sites(user_input, genome)
        # Infer regulons
        prior = get_prior(genome, user_input, weights)
        regulated_genes = infer_regulations(user_input, genome, prior)
        all_regulated_genes.extend(regulated_genes)
        # Output operons
        output_operons(user_input, genome)

    # Create orthologous groups
    grps = create_orthologous_groups(user_input, all_regulated_genes, genomes)
