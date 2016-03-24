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
    """Makes a directory specified by one or more args.

    Creates the directory if it doesn't exist.

    Args:
        paths ([string]): the full path for the directory
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

    Args:
        user_input (UserInput): a UserInput object keeping the user input
    Returns:
        [Genome]: list of created genomes
    """
    my_logger.info("Started: create genomes")
    genomes = [Genome(name, accessions)
               for name, accessions in user_input.genome_name_and_accessions]
    my_logger.info("Finished: create genomes")
    return genomes


def output_operons(user_input, genome):
    """Writes operons of the genome to its corresponding CSV file.

    Args:
        user_input (UserInput): the parameters provided by the user.
        genome (Genome): the genome the operons of which are written to file.
    """
    log_dir = directory(user_input.log_dir, 'operons')
    genome.operons_to_csv(
        os.path.join(log_dir, genome.strain_name+'_operons.csv'))


def identify_TF_instance_in_genomes(genomes, proteins):
    """Identifies TF instance for each genome based on given proteins.

    Args:
        genomes ([Genome]): the list of genomes to identify TF-instance on
        proteins ([Protein]): the list of proteins provided by the user,
            along its binding motifs.
    """
    for g in genomes:
        my_logger.info("Identifying TF instance for /%s/" % g.strain_name)
        g.identify_TF_instance(proteins)


def set_TF_binding_model(user_input, genome, site_collections, weights):
    """Sets the TF-binding model for each genome.

    The TF_binding_models are genome-specific and built using the collections
    of sites and their associated weights.

    For each genome, built models are written to corresponding CSV files.

    Args:
        user_input (UserInput): options provided by the user
        genome (Genome): genome of interest
        site_collections ([SiteCollection]): list of collections
        weights ([float]): weights used for combining evidence
    """
    genome.build_PSSM_model(site_collections, weights)
    log_dir = directory(user_input.log_dir, 'derived_PSWM')
    genome.output_TF_binding_model(
        os.path.join(log_dir, genome.strain_name+'.jaspar'))


def create_proteins(user_input):
    """Creates all proteins used in the analysis.

    Args:
        user_input (UserInput): options provided by the user
    Returns:
        [Protein]: list of created proteins
    """
    my_logger.info("Started: create proteins")
    proteins = [Protein(accession, user_input.TF_name)
                for accession in user_input.protein_accessions]
    my_logger.info("Finished: create proteins")
    return proteins


def create_site_collections(user_input, proteins):
    """Creates SiteCollection objects from the user-provided collections.

    See site_collection.py for SiteCollection definition.

    Args:
        user_input (UserInput): parameters provided by the user
        proteins ([Protein]): list of proteins associated with the collections
    Returns:
        [SiteCollection]: the list of SiteCollection objects, one per protein.
    """
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
    """Weights collections according to collection sizes.

    Args:
        [SiteCollection]: list of SiteCollection objects.
    Returns:
        [float]: list of weights, not normalized
    """
    return [collection.site_count for collection in site_collections]


def phylogenetic_weighting(site_collections, genome, phylogeny):
    """Computes weights of collections based on the phylogeny.

    The weights are inversely proportional to the phylogenetic distances the
    proteins that the motifs belong to, and the TF-instance of the genome of
    interest.

    Args:
        site_collections ([SiteCollection]): list of SiteCollection objects,
            built with the collections of binding sites, provided by the user.
        genome (Genome): the genome of interest
        phylogeny (Phylo): the phlogenetic tree of the instances of the TF of
            interest.
    Returns:
        [float]: list of weights, not normalized.
    """
    dists = [phylogeny.distance(collection.TF, genome.TF_instance)
             for collection in site_collections]
    normalized_dists = [float(d)/sum(dists) for d in dists]
    similarities = [1-nd for nd in normalized_dists]
    weights = [s/sum(similarities) for s in similarities]
    return weights


def uniform_weighting(site_collections):
    """Returns the uniform weights."""
    return [1 for _ in site_collections]


def compute_motif_weights(genome, site_collections, weighting_scheme,
                          phylogeny):
    """Compute motif weights for given genome.

    Two weighting schemes are supported:
    - simple: the weights are computed as proportional to the site collection
      sites.
    - phylogenetic: the weights are genome specific. They are inversely
      proportional to the phylogenetic distance between the genome of interest
      and each protein provided along with the binding evidence.

    Args:
        genome (Genome): the genome the weights for which are computed
        site_collections ([SiteCollection]): list of binding site collections
        weighting_scheme (string): simple or phylogenetic.
        phylogeny (Phylo): the phylogeny to be used to compute the
            weighting. Used only if phylogenetic weighting scheme is set.
    Returns:
        [float]: weights of the site collections
    """
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


def get_prior(genome, user_input, weights):
    """Infers the prior probability of binding.

    This is done by predicting operons in the genome of the TF for which the
    motif is provided and dividing the number of sites in the motif (assuming
    that there is one site per operon) by the number of operons. For instance,
    if 30 sites are available for LexA in E. coli, then the prior for
    regulation is 30/~2300 operons.

    Then it is assigned to genome using the same strategy defined for mixing
    the TF-binding models. That is, simple combination (arithmetic mean of
    priors) or weighted combination (with phylogenetic distances).

    Args:
        genome (Genome): the genome of interest
        user_input (UserInput): the user-provided options.
        weights: weights of each site collection to be compute the final prior
    Returns:
        float: the prior probability of TF binding the DNA.
    """
    # If the prior regulation probability is set by the user, use that value.
    if user_input.prior_regulation_probability:
        return user_input.prior_regulation_probability
    # Otherwise, infer the prior probability.
    my_logger.info("Prior probability not provided, "
                   "inferring from the provided motifs. /%s"
                   % genome.strain_name)
    site_collections = genome.TF_binding_model.site_collections
    priors = [float(collection.site_count)/genome.num_operons
              for collection in site_collections]
    return sum(p*w for p, w in zip(priors, weights))


def infer_regulations(user_input, genome, prior):
    """Scans the given for TF-binding sites and writes results to csv files.

    The operons with a promoter that has TF-binding probability larger than the
    defined threshold (default 0.5) are reported.

    Args:
       user_input (UserInput): the user-provided options
       genome (Genome): the genome of interest
       prior (float): the prior probability of binding to a promoter
    Returns:
       [Gene]: list of regulated genes

    """
    log_dir = directory(user_input.log_dir, 'posterior_probs')
    # Get the probability threshold
    threshold = user_input.probability_threshold
    my_logger.info("Regulation probability threshold: %.2f" % threshold)
    my_logger.info("Inferring regulations for /%s/." % genome.strain_name)
    my_logger.info("Prior probability of regulation: %f" % prior)
    report_filename = os.path.join(log_dir, genome.strain_name+'.csv')
    # Scan the genome to compute binding probability for each promoter.
    genome_scan = genome.infer_regulation(prior, threshold, report_filename)
    # Return the list of regulated genes.
    return [g for opr, _ in genome_scan for g in opr.genes]


def search_sites(user_input, genome):
    """Searches for putative binding sites.

    It uses genome's TF-binding model for binding site search and outputs
    the results into a CSV file accordingly.

    Args:
        user_input (UserInput): the user-provided input
        genome (Genome): genome of interest
    """
    log_dir = directory(user_input.log_dir, 'identified_sites')
    my_logger.info("Scoring genome /%s/." % genome.strain_name)
    my_logger.info(
        "Site score threshold: %.2f" % genome.TF_binding_model.threshold())
    report_filename = os.path.join(log_dir, genome.strain_name+'.csv')
    genome.identify_sites(filename=report_filename)


def create_orthologous_groups(user_input, genes, genomes):
    """Creates orthologous groups from the given list of genes.

    For each gene in the list, it searches for orthologs in other genomes and
    outputs results to the CSV file.

    Args:
        user_input (UserInput): the user-provided input
        genes ([Gene]): the list of genes the orthologs of which are identified
            to create orthologous groups.
        genomes ([Genome]): the list of all genomes to search for orthologs.
    Returns:
        [OrthologousGroup]: the list of orthologous groups.
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

    Args:
        genomes ([Genome]): the list of genomes used in the analysis
        proteins ([Protein]): the list of proteins used in the analysis
        user_input (UserInput): the user-provided input
    Returns:
        Phylo: the phylogeny object. See phylo.py
    """
    phylo = Phylo(proteins + [g.TF_instance for g in genomes])
    # Output the phylogenetic tree in newick format.
    phylo.to_newick(os.path.join(user_input.log_dir, 'phylogeny.nwk'))
    phylo.draw_ascii()
    return phylo


def main():
    """The entry-point for the pipeline."""
    # Read user input and configuration from two files
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
    create_orthologous_groups(user_input, all_regulated_genes, genomes)
