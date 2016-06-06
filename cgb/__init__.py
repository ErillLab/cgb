import os
import pickle
import copy

from .genome import Genome
from .protein import Protein
from .site_collection import SiteCollection
from .my_logger import my_logger
from .misc import directory
from .phylo import Phylo
from .user_input import UserInput
from .orthologous_group import construct_orthologous_groups
from .orthologous_group import orthologous_grps_to_csv
from .orthologous_group import ancestral_state_reconstruction
from .orthologous_group import ancestral_states_to_csv
from .visualization import all_plots


PICKLE_DIR = directory('pickles')
OUTPUT_DIR = directory('output')


def pickle_dump(obj, filename):
    pickle.dump(obj, open(os.path.join(PICKLE_DIR, filename), 'w'))


def create_genomes(user_input):
    """Creates all genomes used in the analysis.

    For each genome, it writes operons to csv files.

    Args:
        user_input (UserInput): a UserInput object keeping the user input
    Returns:
        [Genome]: list of created genomes
    """
    my_logger.info("Started: create genomes")
    # Check for duplicates
    if len(user_input.genome_names) != len(set(user_input.genome_names)):
        my_logger.error("Duplicate in genomes. Exiting...")
        raise SystemExit

    # Create genomes with given names and accession numbers.
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
    output_dir = directory(OUTPUT_DIR, 'operons')
    genome.operons_to_csv(
        os.path.join(output_dir, genome.strain_name+'_operons.csv'))


def identify_TF_instance_in_genomes(genomes, proteins):
    """Identifies TF instance for each genome based on given proteins.

    Args:
        genomes ([Genome]): the list of genomes to identify TF-instance on
        proteins ([Protein]): the list of proteins provided by the user,
            along its binding motifs.
    """
    for g in genomes:
        my_logger.info("Identifying TF instance for (%s)" % g.strain_name)
        g.identify_TF_instance(proteins)


def remove_genomes_with_no_TF_instance(genomes):
    """Removes the genomes with no identified TF-instance from the analysis."""
    genomes_with_TF = [genome for genome in genomes if genome.TF_instance]
    return genomes_with_TF


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
    output_dir = directory(OUTPUT_DIR, 'derived_PSWM')
    genome.output_TF_binding_model(
        os.path.join(output_dir, genome.strain_name+'.jaspar'))


def create_proteins(user_input):
    """Creates all proteins used in the analysis.

    Args:
        user_input (UserInput): options provided by the user
    Returns:
        [Protein]: list of created proteins
    """
    my_logger.info("Started: create proteins")
    proteins = []
    for accession in user_input.protein_accessions:
        my_logger.info("Initializing %s." % accession)
        proteins.append(Protein(accession, user_input.TF_name))
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
    output_dir = directory(OUTPUT_DIR, 'user_PSWM')
    for collection in collections:
        collection.to_jaspar(
            os.path.join(output_dir, collection.TF.accession_number+'.jaspar'))
    my_logger.info("Finished: create site collections")
    return collections


def site_count_weighting(site_collections):
    """Weights collections according to collection sizes.

    Args:
        [SiteCollection]: list of SiteCollection objects.
    Returns:
        [float]: list of weights, not normalized
    """
    return [collection.site_count for collection in site_collections]


def phylogenetic_weighting(site_collections, genome, phylogeny,
                           clustalesque_weighting=True):
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
        clustalesque_weighting (bool): if true, the weights are still computed
            based on the phylogenetic distance, but branch lengths are inflated
            proportional to the number of terminal nodes in each branch. Doing
            so would down-weight the evidence from closer species and up-weight
            the evidence from the most divergent ones.
    Returns:
        [float]: list of weights, not normalized.
    """

    # Collect the protein accession numbers with binding evidence.
    reference_proteins = [collection.TF.accession_number
                          for collection in site_collections]
    # Don't modify the original tree
    tree = copy.deepcopy(phylogeny.tree)
    # Prune nodes except the target species and the species with evidence
    for c in tree.get_terminals():
        if not (c.name == genome.TF_instance.accession_number or
                c.name in reference_proteins):
            tree.prune(c)

    # Convert phylogenetic distances to similarities
    total_branch_length = tree.total_branch_length()
    for c in tree.find_clades():
        if c.branch_length:
            c.branch_length = 1.0 - c.branch_length / total_branch_length

    if clustalesque_weighting:
        # Weight similarities like CLUSTAL does for multiple sequence alignment
        # Root the tree using the reference TF as outgroup
        outgroup = tree.find_any(name=genome.TF_instance.accession_number)
        tree.root_with_outgroup(outgroup)
        # Down-weight all branches with multiple species with evidence
        for protein in reference_proteins:
            node = tree.find_any(name=protein)
            for c in tree.get_path(node):
                c.branch_length /= c.count_terminals()

    node = tree.find_any(name=genome.TF_instance.accession_number)
    weights = []
    for acc in reference_proteins:
        other = tree.find_any(name=acc)
        weights.append(tree.distance(node, other))

    return weights


def compute_weights(genome, site_collections, phylogeny=None,
                    site_counts=False):
    """Computes weights of each binding evidence for the target genome.

    Args:
        genome (Genome): the target genome that the weights are for
        site_collections ([SiteCollection]): list of binding site collections
        phylogeny (Phylo): the phylogeny to be used for phylogenetic
            weighting. If none, no phylogenetic weighting is applied.
        site_counts (bool): if true, evidence is weighted by the number of
            sites in each collection
    """
    my_logger.info("Computing weights for %s..." % genome.strain_name)
    my_logger.info("Use phylogeny: %s" % bool(phylogeny))
    my_logger.info("Use site counts: %s" % bool(site_counts))

    my_logger.info('Site collections from %s' %
                   [c.TF.accession_number for c in site_collections])
    weights = [1.0 for _ in site_collections]  # Equal weights by default
    if phylogeny:
        phylo_ws = phylogenetic_weighting(site_collections, genome, phylogeny)
        my_logger.info('Phylogeny weights: %s' %
                       ['%.2f' % w for w in phylo_ws])
        # Update weights
        weights = [w*pw for w, pw in zip(weights, phylo_ws)]
    if site_counts:
        site_count_ws = site_count_weighting(site_collections)
        my_logger.info('Site count weights: %s' %
                       ['%.2f' % w for w in site_count_ws])
        # Update weights
        weights = [w*cw for w, cw in zip(weights, site_count_ws)]

    # Normalize weights
    weights = [w/sum(weights) for w in weights]
    my_logger.info('Final weights: %s' % ['%.2f' % w for w in weights])

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
                   "inferring from the provided motifs. (%s)"
                   % genome.strain_name)

    print genome.TF_binding_model.IC
    prior = (genome.length /
             2**genome.TF_binding_model.IC /
             genome.num_operons)
    my_logger.info("Prior for %s: %f" % (genome.strain_name, prior))
    return prior


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
    output_dir = directory(OUTPUT_DIR, 'posterior_probs')
    # Get the probability threshold
    threshold = user_input.probability_threshold
    my_logger.info("Regulation probability threshold: %.2f" % threshold)
    my_logger.info("Inferring regulations for (%s)." % genome.strain_name)
    my_logger.info("Prior probability of regulation: %f" % prior)
    report_filename = os.path.join(output_dir, genome.strain_name+'.csv')
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
    output_dir = directory(OUTPUT_DIR, 'identified_sites')
    my_logger.info("Scoring genome (%s)." % genome.strain_name)
    my_logger.info(
        "Site score threshold: %.2f" % genome.TF_binding_model.threshold())
    report_filename = os.path.join(output_dir, genome.strain_name+'.csv')
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
    orthologous_grps_to_csv(groups, os.path.join(OUTPUT_DIR, 'orthologs.csv'))
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
    print phylo.tree
    # Output the phylogenetic tree in newick and nexus formats.
    phylo.to_newick(os.path.join(OUTPUT_DIR, 'phylogeny.nwk'))
    phylo.to_nexus(os.path.join(OUTPUT_DIR, 'phylogeny.nex'))
    return phylo


def perform_ancestral_state_reconstruction(user_input, genomes,
                                           orthologous_grps):
    """Performs ancestral state reconstruction

    It uses BayesTraitsV2 (http://www.evolution.rdg.ac.uk/BayesTraits.html) to
    estimate trait probabilities on intermediate nodes of the given
    phylogenetic tree.

    It outputs the estimated state probabilities to a CSV file.

    It draws the phylogenetic tree of target genomes to a file.

    Args:
        user_input (UserInput): the parameters provided by the user
        genomes ([Genome]): the list of target genomes
        orthologous_grps ([OrthologousGroup]): the list of orthologous groups
            for which the ancestral state reconstruction is performed. Each
            group consists of genes that are orthologous to each other.
    """
    # Create a phylogeny using target genomes only.
    phylo = Phylo([g.TF_instance for g in genomes],
                  names=[g.strain_name for g in genomes])
    # Perform ancestral state reconstruction
    ancestral_state_reconstruction(orthologous_grps, phylo)

    # Write ancestral state reconstruction to a csv file
    ancestral_states_to_csv(orthologous_grps, phylo,
                            os.path.join(OUTPUT_DIR, 'ancestral_states.csv'))
    # Save phylogeny
    phylo.draw(os.path.join(OUTPUT_DIR, "phylogeny.png"))


def go(input_file):
    """The entry-point for the pipeline."""
    # Read user input and configuration from two files
    user_input = UserInput(input_file)
    pickle_dump(user_input, 'user_input.pkl')

    # Make output directory
    directory(OUTPUT_DIR)

    # Create proteins
    proteins = create_proteins(user_input)
    pickle_dump(proteins, 'proteins.pkl')

    # Create genomes
    genomes = create_genomes(user_input)

    # Identify TF instances for each genome.
    identify_TF_instance_in_genomes(genomes, proteins)
    # Remove genomes with no TF-instance from the analysis.
    genomes = remove_genomes_with_no_TF_instance(genomes)

    # Create phylogeny
    phylogeny = create_phylogeny(genomes, proteins, user_input)
    pickle_dump(phylogeny, 'phylogeny.pkl')

    # Create binding evidence
    site_collections = create_site_collections(user_input, proteins)

    all_regulated_genes = []
    for genome in genomes:
        # Create weights for combining motifs and inferring priors, if not set
        phylo_weighting = (phylogeny if user_input.phylogenetic_weighting
                           else None)
        site_count_weighting = user_input.site_count_weighting
        weights = compute_weights(genome, site_collections,
                                  phylo_weighting, site_count_weighting)
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
    pickle_dump(genomes, 'genomes.pkl')

    # Create orthologous groups
    ortholog_groups = create_orthologous_groups(
        user_input, all_regulated_genes, genomes)
    pickle_dump(ortholog_groups, 'orthos.pkl')

    # Ancestral state reconstruction step
    perform_ancestral_state_reconstruction(
        user_input, genomes, ortholog_groups)

    pickle_dump(ortholog_groups, 'orthos2.pkl')

    # Create phylogenetic tree of target genomes only
    phylo_target_genomes = Phylo([g.TF_instance for g in genomes],
                                 names=[g.strain_name for g in genomes])
    all_plots(phylo_target_genomes, ortholog_groups, genomes,
              directory(OUTPUT_DIR, 'plots'))
