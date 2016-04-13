import subprocess
import os
import re

from misc import temp_file_name

# Path to BayesTraits executable
BAYES_TRAITS = './BayesTraitsV2'


def generate_tree_file(phylo):
    """Generates the tree file that BayesTraits uses.
    Args:
        phylo (Phylo): the phylogeny object
    Returns:
        string: the Nexus file that the tree is written into.
    """
    tree_file = temp_file_name()
    phylo.to_nexus(tree_file)
    return tree_file


def generate_traits_file(trait):
    """Generates the data file containing traits data.

    Args:
        trait (dict): dictionary in {species name: state} format.
    Returns:
        string: the data file to be used by BayesTraits.
    """
    data_file = temp_file_name()
    with open(data_file, 'w') as f:
        for sp, state in trait.items():
            f.write('%s %s\n' % (sp, state))
    return data_file


def generate_command_file(phylo):
    """Generates the command file to be passed to the BayesTraits.

    Args:
        phylo (Phylo): the phylogeny object
    Returns:
        string: the command file to be passed to BayesTraits.
    """
    command_file = temp_file_name()
    with open(command_file, 'w') as f:
        f.write("1\n")              # Select MultiState
        f.write("1\n")              # Select Maximum-Likelihood
        # Create internal nodes for BayesTraits
        for node in phylo.tree.get_nonterminals():
            f.write("AddNode %s %s\n" % (
                node.name, ' '.join(n.name for n in node.get_terminals())))

        f.write("Run\n")            # Run command
    return command_file


def run_bayes_traits(tree_file, trait_file, command_file):
    """Runs BayesTraits.

    Calls BayesTraits with the given tree file, data file and command file.
    The results are stored in a file named as {trait_file}.log.txt

    Args:
        tree_file (string): name of the tree file
        trait_file (string): name of the trait data file
        command_file (string): name of the command file
    """
    subprocess.call([BAYES_TRAITS, tree_file, trait_file],
                    stdin=open(command_file),
                    stdout=open(os.devnull, 'w'))


def parse_bayes_trait_output(result_file):
    """Parses BayesTraits output file.

    Args:
        result_file (string): the name of the result file to be parsed
    Return:
        dictionary containing probability of each state for each ancestral
        node. The dictionary has {(node name, state): probability} format.
    """
    # Read result file
    with open(result_file) as f:
        lines = f.read().splitlines()

    if lines[0].startswith("There has to be more then one state in file"):
        # All species have the same value for the trait. BayesTraits doesn't
        # perform ancestral state reconstruction
        raise OneStateException

    results = {}
    for field, val in zip(lines[-2].split('\t'), lines[-1].split('\t')):
        match = re.match(r'(\S+) P\((\w)\)', field)
        if match:
            node, state = match.group(1), match.group(2)
            results[(node, state)] = float(val)
    return results


def bayes_traits(phylo, trait):
    """Entry point for BayesTraits.

    Runs the BayesTraitsV2 with the given phylogenetic tree and states on
    leaves, parses estimated ancestral state probabilities from the BayesTraits
    output file and returns the ancestral states.

    Args:
        phylo (Phylo): the phylogeny of target genomes
        trait ({species: state}): the dictionary containing the regulation
            state for each species.
    Returns:
        A dictionary containing probability of each state for each ancestral
        node. The dictionary has {(node_name, state): probability} format.
    """
    tree_file = generate_tree_file(phylo)
    data_file = generate_traits_file(trait)
    command_file = generate_command_file(phylo)

    run_bayes_traits(tree_file, data_file, command_file)
    # Check whether all traits have the same state. If so, all ancestral states
    # would have the same trait with probability 1.0
    try:
        log_file = data_file + '.log.txt'
        ancestral_states = parse_bayes_trait_output(log_file)
    except OneStateException:
        # It means all species have the same value for the trait. In this case,
        # BayesTraits doesn't perform ancestral state reconstruction as it
        # requires at least two states in the data file. Here, it is assumed
        # that all ancestral states would have the only state with probability
        # of 1.
        state, = set(trait.values())
        ancestral_states = {(node.name, state): 1.0
                            for node in phylo.tree.get_nonterminals()}
    return ancestral_states


class OneStateException(Exception):
    """Raised when all species have the same value for the trait."""
    pass
