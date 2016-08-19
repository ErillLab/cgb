CGB outline
===========

The overall structure of CGB is a fairly linear pipeline: reading input files and preprocessing, performing actual search & applying comparative genomics filtering to them.

Input reading 
-------------

Input is entered via a **JSON file**, sitting on the root folder of the repo, which defines the following (as lists of dictionaries):

* motifs: a number of experimentally-determined TF-binding sites, together with the protein accession for the TF and the spcecies name (used for the tree)
* genomes: all the genomes to be analyzed, composed of a name and a bunch of accession numbers (can be draft/incomplete genomes)
* parameters:

  * prior_regulation_probability: genomic prior for *operon* regulation
  * phylogenetic_weighting: ON/OFF
  * site_count_weighting: ON/OFF
  * posterior_probability_threshold: posterior cut for regulation assignment
  * posterior_probability_threshold_for_operon_prediction: posterior cut for splitting operons
  * ancestral_state_reconstruction: ON/OFF
  * operon_prediction_distance_tuning_parameter: used to increase/decrease the computed mean intergenic distance
  

The JSON file input is saved onto a **UserInput** Python object, which essentially provides get methods.

After reading the JSON input, the pipeline calls **create_proteins**, which creates protein objects containing the NCBI protein record, and **create_genomes** which does the same thing for genomes. The are stored in the *proteins* and *genomes* objects.

Preprocessing
-------------

TF detection
____________

Next, the pipeline invokes **identify_TF_instance_in_genomes**, which will use BLAST to identify instances of each TF in each of the genomes for analysis. To do this, the system uses the **genome** class method **identify_TF_instance**, which takes protein objects (all the TF instances provided) and tries to identify them through BLAST (with eval=10**-3) in each genome. If there are any genomes with no orthologs of the TF, they will be dropped from the analysis through **remove_genomes_with_no_TF_instance**, which essentially removes them from the *genomes* list.

Phylogeny
_________

A phylogeny including the mapped TFs and the reference ones is created by **create_phylogeny(genomes, proteins, user_input)**, which uses methods from the CGB phylo.py module. The phylo.py module is a wrapper around the BioPython Phylo and Align libs.

The main function in phylo.py is **tree**, which computes the phylogenetic tree. To do so, it computes all the pair-wise protein distances (as %identities), invokes the NJ method from the BioPhylo library and roots it at mid-poiint.

Site collections
________________

Next, the pipeline creates the motifs (i.e. site collections) for the TF instances provided by the user, using **create_site_collections**. Site collections are stored in the **collections** object. Each site collection has an associated TF, the site instances, and a BioPython motifs object.

Regulon prediction 
------------------
The main part of the CGB pipeline is devoted to identifying putative TF-binding sites in all genomes, thus inferring a putative regulon for the TF of interest. This s done in a loop (one genome at a time) and involves some complex operations, having to do mainly with two things: how to define the TF-binding model in each genome (i.e. combining experimental models) and how to define operons.

Weight computation
__________________

The first step in the main pipeline loop is to compute the weighting scheme to transfer information into the genome. This information is, essentially, the TF-binding model and the prior for regulation.

The function that computes the weights is **compute_weights**. It takes in the genome, all site collections, the phylogeny (or none if not used) and a boolean indicator for site count weighting. For phylogeny weighting, the function calls **phylogenetic_weighting** and returns a set of weights for each site collection. The same process is used for collection size, using **site_count_weighting**. If no weighting scheme is used, the weights will be one by default. The weights are multiplicative factors on the additive normalized sum of the contributions to the TF-model or prior probabilty for that genome.

The **phylogenetic_weighting** function incorporates a *"clustalesque"* weighting scheme. The idea here is the same as in CLUSTALW. If a branch of the tree has a lot of evidence, and another only a bit, the "populated" branch is going to dominate, and set all motifs and priors according to that side of the tree (even in spite of contradictory evidence somewhere else on the tree). To counteract this, the *clustalesque* approach dillutes the contribution of a populated branch by recomputing the distance of the elements in that branch to the target genome, essentially making them "more distant" (proportionally to how many informative elements lie below the branch).

In the case of **phylogenetic_weighting**, the tree is first rerooted using the genome of interest as the outgroup. Then, the  branch length of each clade (what connects it to the upper part of the tree) is multiplied by the number of informative elements within the clade (this pushes the clade away proportionally to how much information-bearing species it contains). This is done for all clades.

After this is done (if option enabled), weighting ensues. The target genome is identified on the tree, and for each site collection, we find its node on the tree and measure the distance between the target genome and the site collection node. These distances are then transformed to weights, by computing [1.0 - (dist - min(distances))/max(distances). Weights are then normalized to add up to one. The (1-min)/max seeks to get a 1 (maximal) contribution if the distance between target and refernece is small (min), or a very small contribution (max-min)/max if the reference is very far away (max).
