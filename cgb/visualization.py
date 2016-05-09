
import os
import copy
import math

from ete3 import Tree
from ete3 import TreeStyle, NodeStyle
from ete3 import RectFace, TextFace, CircleFace, StackedBarFace, ImgFace
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from my_logger import my_logger


def biopython_to_ete3(biopython_tree):
    """Converts Biopython Bio.Phylo.Tree to ETE3 TreeNode.

    ETE allows customizing the visualization of phylogenetic
    trees. Documentation for ETE toolkit is available at
    http://etetoolkit.org/docs/latest/index.html
    """
    t = Tree(biopython_tree.format('newick'), format=1)
    return t


def rgb2hex(red, green, blue):
    """Converts the given rgb values to hexadecimal color code.

    Assumes that rgb values are positive and sum to 1.
    """
    return '#%02x%02x%02x' % (red*255, green*255, blue*255)


def filter_and_sort_orthologous_grps(orthologous_groups, min_size=2):
    orthos = [grp for grp in orthologous_groups if len(grp.genes) >= min_size]

    orthos.sort(
        key=lambda grp:
        len([gene for gene in grp.genes
             if grp.most_likely_state_at(gene.genome.strain_name) == '1']),
        reverse=True)
    return orthos


def heatmap_view(tree, orthologous_groups, save_dir):
    """Generates a heatmap of regulation states in all species."""
    # Heat map settings
    rect_face_fgcolor = 'black'
    locus_tag_len = max(len(gene.locus_tag)
                        for ortho_grp in orthologous_groups
                        for gene in ortho_grp.genes)
    rect_face_width = locus_tag_len * 8
    rect_face_height = 20
    rotation = 90

    # Sort orthologous groups by the number of regulated genes in each group
    orthologous_groups = filter_and_sort_orthologous_grps(orthologous_groups)

    # For each species and its gene in each orthologous group, draw a rectangle
    for node in tree.get_leaves():
        for i, orthologous_grp in enumerate(orthologous_groups, start=1):
            reg_states = orthologous_grp.regulation_states
            p_regulation = reg_states[(node.name, '1')]
            p_notregulation = reg_states[(node.name, '0')]
            p_absence = reg_states[(node.name, 'A')]
            # Color of the rectangle is based on probabilities
            rect_face_bgcolor = rgb2hex(p_notregulation, p_regulation, p_absence)
            # Check if this species have a gene in this orthologous group
            gene = orthologous_grp.member_from_genome(node.name)
            rect_face_text = gene.locus_tag if gene else ''
            rect_face_label = {'text': rect_face_text,
                               'font': 'Courier',
                               'fontsize': 8,
                               'color': 'black'}
            # Create the rectangle
            rect_face = RectFace(rect_face_width, rect_face_height,
                                 rect_face_fgcolor, rect_face_bgcolor,
                                 label=rect_face_label)

            rect_face.rotation = -rotation
            # Add the rectangle to the corresponding column
            node.add_face(rect_face, column=i, position='aligned')

    ts = TreeStyle()
    # Add orthologous group descriptions
    descriptions = [grp.description for grp in orthologous_groups]
    max_description_len = max(map(len, descriptions))
    descriptions = [
        '(%d) ' % i + description + ' '*(max_description_len-len(description))
        for i, description in enumerate(descriptions, start=1)]
    # Append space to make descriptions look nicer

    for i, description in enumerate(descriptions, start=1):
        text_face = TextFace(description, ftype='Courier')
        text_face.hz_align = 1
        text_face.vt_align = 1
        text_face.rotation = -rotation
        ts.aligned_header.add_face(text_face, column=i)

    # Rotate the generated heatmap.
    ts.margin_left = 10
    ts.margin_top = 20
    ts.rotation = rotation
    # For some reason, it can't render to PDF in color
    tree.render(os.path.join(save_dir, 'heatmap.svg'), tree_style=ts)


def view_by_gene(tree, orthologous_grp, save_file):
    # Make a copy of the tree
    tree = copy.deepcopy(tree)
    for node in tree.traverse():
        percents = [orthologous_grp.regulation_states[(node.name, s)] * 100
                    for s in ['1', '0', 'A']]
        stacked_bar_face = StackedBarFace(percents, width=30, height=15,
                                          colors=['green', 'red', 'blue'])
        stacked_bar_face.margin_left = 5
        stacked_bar_face.margin_right = 5
        node.add_face(stacked_bar_face, column=1)
         # Add locus tag if the genome has the gene from this orthologous group
        if node.is_leaf():
            # Check if the orthologous group contains any gene of the genome
            gene = orthologous_grp.member_from_genome(node.name)
            if gene:
                node.add_face(TextFace(text='[%s]' % gene.locus_tag), column=2)

    ts = TreeStyle()
    ts.title.add_face(TextFace(text=orthologous_grp.description), column=1)
    tree.render(save_file, tree_style=ts)


def view_all_genes(tree, orthologous_groups, save_dir):
    orthologous_groups = filter_and_sort_orthologous_grps(orthologous_groups)
    for i in tqdm(range(len(orthologous_groups))):
        view_by_gene(tree, orthologous_groups[i],
                     os.path.join(save_dir, 'ortho_%04d.svg' % i))


def view_by_taxon_(tree, orthologous_groups, node_name, display_threshold=0.95):
    graph = nx.DiGraph()
    graph.add_node(0, color='yellow', label='TF')

    for i, orthologous_grp in enumerate(orthologous_groups, 1):
        p_reg = orthologous_grp.regulation_states[(node_name, '1')]
        p_not_reg = orthologous_grp.regulation_states[(node_name, '0')]
        p_absence = orthologous_grp.regulation_states[(node_name, 'A')]
        if p_reg <= display_threshold:
            continue
        description = (orthologous_grp.member_from_genome(node_name).product
                       if orthologous_grp.member_from_genome(node_name)
                       else orthologous_grp.description)
        graph.add_node(i, fillcolor=rgb2hex(p_not_reg, p_reg, p_absence),
                       label=description)
        graph.add_edge(0, i, weight=p_reg)

    graph.add_edge(0, 0, weight=1)

    #print graph.nodes()
    #pos = nx.drawing.nx_agraph.graphviz_layout(graph)
    #nx.draw_networkx_nodes(graph, pos=pos, node_color=colors, alpha=0.8)
    #nx.draw_networkx_edges(graph, pos=pos, alpha=0.2)
    #nx.draw_networkx_labels(graph, pos=pos)
    #nx.draw_networkx_labels(graph, pos=pos, labels=labels)

    agraph = nx.nx_agraph.to_agraph(graph)
    agraph.node_attr['style'] = 'filled'
    agraph.node_attr['shape'] = 'rectangle'
    agraph.draw('file.png', prog='circo')


def view_by_taxon(node_name, orthologous_groups, save_file):
    # Find the orthologous groups that contains or are likely to contain a gene
    # of the taxon.
    regulation_probs = []
    not_regulation_probs = []
    absence_probs = []
    locus_tags = []
    products = []
    # Sort orthologous groups by regulation probability at the taxon
    orthologous_groups.sort(key=lambda grp: grp.regulation_states[(node_name, '1')])
    for ortho_grp in orthologous_groups:
        reg_prob = ortho_grp.regulation_states[(node_name, '1')]
        not_reg_prob = ortho_grp.regulation_states[(node_name, '0')]
        absence_prob = ortho_grp.regulation_states[(node_name, 'A')]
        if absence_prob > 0.99:
            continue
        regulation_probs.append(reg_prob)
        not_regulation_probs.append(not_reg_prob)
        absence_probs.append(absence_prob)
        gene = ortho_grp.member_from_genome(node_name)
        if gene:
            locus_tags.append(gene.locus_tag)
            products.append(gene.product)
        else:
            locus_tags.append('')
            products.append(ortho_grp.description)

    height = 0.5
    N = len(regulation_probs)

    ind = np.arange(len(regulation_probs))
    fig = plt.figure(figsize=(10, 0.6*N))
    ax1 = fig.add_subplot(111)
    rects1 = ax1.barh(ind, regulation_probs, height=height, color='g')
    rects2 = ax1.barh(ind, not_regulation_probs, height=height, color='r',
                      left=regulation_probs)
    rects3 = ax1.barh(
        ind, absence_probs, color='b', height=height,
        left=map(sum, zip(regulation_probs, not_regulation_probs)))
    ax1.set_xlabel('Probability')
    ax1.set_yticks(ind+height/2.)
    ax1.set_yticklabels(locus_tags)
    ax1.set_ylim(0, N+2)

    ax2 = ax1.twinx()
    ax2.set_yticks(ax1.get_yticks())
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_yticklabels(products)

    ax1.legend((rects1[0], rects2[0], rects3[0]),
               ('Regulation', 'Not regulation', 'Absence'),
               shadow=True, loc='best')

    plt.savefig(save_file, bbox_inches='tight')


def view_all_taxa(tree, orthologous_grps, save_dir):
    for node in tree.traverse():
        my_logger.info("Plotting %s." % node.name)
        view_by_taxon(node.name, orthologous_grps,
                      os.path.join(save_dir, 'taxon_%s.svg' % node.name))


def network_size_view(tree, orthologous_groups, save_dir):
    # Set the size of each node
    for node in tree.traverse():
        # Count the number of genes likely to be regulated in this node
        num_regulated = len([grp for grp in orthologous_groups
                             if grp.most_likely_state_at(node.name) == '1'])
        nstyle = NodeStyle()
        nstyle['size'] = 0
        node.set_style(nstyle)
        circle_face = CircleFace(
            radius=num_regulated, color='green', style='sphere',
            label={'text': str(num_regulated), 'color': 'black'})
        node.add_face(circle_face, column=2)
    ts = TreeStyle()
    tree.render(os.path.join(save_dir, 'network_size.svg'), tree_style=ts)


def motif_view(tree, genomes, save_dir):
    for node in tree.get_leaves():
        genome, = [genome for genome in genomes
                   if genome.strain_name == node.name]
        motif_file = genome.TF_binding_model.weblogo_from_pwm
        img_face = ImgFace(motif_file)
        node.add_face(img_face, column=1, position='aligned')
    ts = TreeStyle()
    tree.render(os.path.join(save_dir, 'binding_motif.svg'), tree_style=ts)


def all_plots(phylo, orthologous_groups, genomes, save_dir):
    """Draw all plots and save them in save_dir."""
    # Heat map
    heatmap_view(biopython_to_ete3(phylo.tree), orthologous_groups, save_dir)

    # Orthologous group-centric view
    view_all_genes(biopython_to_ete3(phylo.tree), orthologous_groups, save_dir)

    # Taxa-centric view
    view_all_taxa(biopython_to_ete3(phylo.tree), orthologous_groups, save_dir)

    # Network size plot
    network_size_view(biopython_to_ete3(phylo.tree), orthologous_groups, save_dir)

    # Motif phylogeny plot
    motif_view(biopython_to_ete3(phylo.tree), genomes, save_dir)
