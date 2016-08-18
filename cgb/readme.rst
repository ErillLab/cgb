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
The main part of the CGB pipeline is devoted to identifying putative TF-binding sites in all genomes, thus inferring a putative regulon for the TF of interest. This involves some complex operations, having to do mainly with two things: how to define the TF-binding model in each genome (i.e. combining experimental models) and how to define operons.

