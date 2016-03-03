import json
import os

from genome import Genome
from protein import Protein
from site_collection import SiteCollection
from my_logger import my_logger
from phylo import Phylo

def parse_input(filename):
    """Parses the input file and returns a parameter dictionary."""
    my_logger.info("Reading input from %s." % filename)
    with open(filename) as f:
        data = json.load(f)
    return data

def create_genomes(genome_input):
    my_logger.info("Started: create genomes")
    genomes = [Genome(g['name'], g['accession_numbers']) for g in genome_input]
    my_logger.info("Finished: create genomes")
    return genomes


def create_proteins(input):
    my_logger.info("Started: create proteins")
    proteins = [Protein(m['protein_accession'], input['TF'])
                for m in input['motifs']]
    my_logger.info("Finished: create proteins")
    return proteins


def get_site_lists(motif_input):
    return [m['sites'] for m in motif_input]


def create_site_collections(proteins, input):
    my_logger.info("Started: create site collections")
    collections = get_site_lists(input['motifs'])
    site_collections = [SiteCollection(sites, protein)
                        for sites, protein in zip(collections, proteins)]
    my_logger.debug("Writing site collections.")
    log_dir = os.path.join(input['configuration']['log_dir'], 'user_PSWM')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    for collection in site_collections:
        collection.to_jaspar(
            os.path.join(log_dir, collection.TF.accession_number+'.jaspar'))
    my_logger.info("Finished: create site collections")
    return site_collections


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


def main():
    # Input collection
    input = parse_input('../tests/input.json')
    # Create proteins
    proteins = create_proteins(input)
    # Create genomes
    genomes = create_genomes(input['genomes'])
    # Create binding evidence
    site_collections = create_site_collections(proteins, input)
    collection_weights = compute_motif_weights(site_collections, 'simple')
    for g in genomes:
        # Identify TF instances.
        g.TF_instance, e_val = g.identify_TF_instance(proteins)
        # Set TF-binding model for each genome.
        prior_reg = 0.01
        g.build_PSSM_model(site_collections, collection_weights, prior_reg)
        log_dir = os.path.join(input['configuration']['log_dir'], 'derived_PSWM')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        g.PSSM_model_to_jaspar(os.path.join(log_dir, g.strain_name+'.jaspar'))
        # Compute operons
        log_dir = os.path.join(input['configuration']['log_dir'], 'operons')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        g.operons_to_csv(os.path.join(log_dir, g.strain_name+'_operons.csv'))

    phylo = Phylo(proteins + [g.TF_instance for g in genomes])
    phylo.to_newick(os.path.join(input['configuration']['log_dir'], 'phylogeny.nwk'))
