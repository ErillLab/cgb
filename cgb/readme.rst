CGB outline
===========

The overall structure of CGB is a fairly linear pipeline: reading input files and preprocessing, performing actual search & applying comparative genomics filtering to them.

Input reading and preprocessing
-------------------------------

Input is entered via a JSON file, sitting on the root folder of the repo, which defines the following (as lists of dictionaries):

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
  


