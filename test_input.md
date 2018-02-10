# CGB Input File overview

CGB input is specified in JSON format, which essentially defines elements as dictionaries and lists thereof. In the case of CGB, the JSON file contains three primary entries (TF, motifs and genomes), and a flurry of configuration parameters. The following sections detail each of these.

## TF
The `TF` field simply specifies the transcription factor for which the analysis is being performed.

`"TF" : "transcription_factor_name_here"`

## Motifs
The `motifs` field defines a list of motifs, each of them defined as a dictionary structure.

`"motifs" : [{...}, {....}]`

### Motif
Each motif is defined as a dictionary that contains the following items:

- Name
	- The name given to the motif. By convention, `TFname_SpeciesAbbr`
	- `"name" : "TF_SpecAbb"`
- Genome accession
	- The list of NCBI accession numbers for the genome (including version number)
	- `"genome_accessions": ["YY_XXXXXX.Z","MM_NNNNNN.O"]`
- Sites
	- A list of aligned transcription factor-binding sites reported (or inferred by the user) for the TF under analysis reported in the genome accessions specified for this organism
	- `"sites" : ["ATGCATCG", "ATCGGCTA"]`

## Genomes
The genomes field defines a list that contains the target genomes for the comparative genomics analysis. Each genome is defined by a dictionary structure.

`"genomes" : [{...}, {...}]`

### Genome
Each genome is defined as a dictionary that contains the following items:

- Name
	- The name given to the target species
	- `"name" : "species_name"`
- Accession numbers
	- A list of NCBI accession numbers for the genome (including version number)
	- `"accession_numbers": ["YY_XXXXXX.Z","MM_NNNNNN.O"]`

## Configuration parameters

### Site search and operon prediction methods
Several parameters control differenet aspects of CGB behavior regarding how it searches for sites and predicts operons.

- `prior_regulation_probability`
	- Type: `mixed`
	- Allowed values:
		- Positive real number in [0..1] interval
		- "No" string
    - Default value: `"No"`
	- Recommended value: `0.03`
		- Alternatively, any estimate by the user on a reference species, using the number of experimentally reported binding sites as the putative number of regulated operons, and dividing this by the number of known (or predicted) operons in that reference genome.
	- Rationale:
		- This parameter allows the user to define a the prior probability of regulation, which is used for the computation of the posterior probability of regulation in CGB.
		- When "No" is specified, CGB will estimate the prior probability of regulation independently in each genome, based on the information content of the TF-binding motif transferred to that genome as a proxy for the predicted number of regulated operons.

- `alpha`
	- Type: `float`
	- Allowed values:
		- Positive real number in [0..1] interval
    - Default value: `0.00285714285714` [1/300]
    - Rationale:
    	- CGB defines the assumes that, in the promoter region of a regulated gene, the distribution of PSSM scores is a mixture of two distributions: the background score distribution and the distribution of scores in functional sites. The former can be inferred from the genome and the latter for the collection (or more generally, PSWM model) of known (or transferred) sites.
    	- For the mixture distribution `R`, the `alpha` parameter defines the mixing ratio between model (`M`) and background (`B`) distributions, so that `R ~ alpha*M + (1-alpha)*B`.
    	- For a *regular* bacterial transcription factor, we expect 1 site per regulated promoter, and an average promoter length of 300 bp, so the default prior is set to 1/300.

- `promoter_up_distance`
	- Type: `integer`
	- Allowed values:
		- Positive integer under 1,000
    - Default value: `300`
    - Rationale:
    	- This parameter controls how far upstream of a gene's translational start site (TLS) CGB will scan for putative TF-binding sites.

- `promoter_dw_distance`
	- Type: `integer`
	- Allowed values:
		- Positive integer under 1,000
    - Default value: `50`
    - Rationale:
    	- This parameter controls how far downstream of a gene's  translational start site (TLS)  CGB will scan for putative TF-binding sites.

- `phylogenetic_weighting`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Recommended value: `true`
    - Rationale:
    	- This paramater determines whether phylogenetic weighting is activated. If activated, the genetic distance between species (estimated from a CLUSTALO) alignment of the homologous TF protein sequences in each target and reference species) is used to weigh how much each reference TF-binding motif contributes to the inferred TF-binding model in each target species.

- `site_count_weighting`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- This paramater determines whether site-count weighting is activated. If activated, the number of sites reported in each reference species is used to weigh the contribution to the inferred TF-binding model in target genomes of different reference species (weighing them up contributions from species with many sites, and down for species with few sites). The aim of this weighting is to counteract artifactual changes in the TF-binding motif due to undersampling.

- `posterior_probability_threshold`
 	- Type: `float`
	- Allowed values:
		- Positive real number in [0..1] interval
    - Default value: `0.5`
    - Rationale:
    	- This parameter defines the minimum posterior probability value required to split an operon. Predicted operons can be split if one of their "internal" genes appears to be regulated. This parameter determines how high the posterior probability of regulation must be in order to split the operon.

### Ancestral state reconstruction

Two parameters control the ancestral state reconstruction process.

- `ancestral_state_reconstruction`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- This paramater determines whether ancestral state reconstruction will be performed.

- `bootstrap_replicates`
	- Type: `integer`
	- Allowed values:
		- Positive integer under 10,000
    - Default value: `100`
    - Recommended value: `100` for fast, preliminary analysis; `1,000` for robust ancestral state reconstruction
    - Rationale:
    	- This parameter sets the number of boostrap replicates generated when performing ancestral state reconstruction. Since ancestral state reconstruction operates with discrete traits, for any given gene the inferred probability of regulation in each species is discretized into 1/0 values (or A for absence of ortholog). This is done `bootstrap_replicates` times, and ancestral reconstruction is run on each replicate. The resulting values for ancestral states in each non-terminal node are averaged.

### Plots and output

Several parameters control the plots and output files generated by CGB.

- `heatmap_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Generates a heatmap plot, with a species phylogram on the x axis and the ortholog groups on the y axis. The color scheme indicates the average posterior probability of regulation for the gene in each species.

- `motif_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- Plots the inferred TF-binding model in each target species using WebLogo.

- `gene_regulation_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Plots the reconstructed ancestral states as a phylogram for each orthologous group.

- `taxon_regulation_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Plots the reconstructed ancestral states for each taxon.

- `network_size_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Plots a phylogram depicting the size of the network in each taxon.

- `site_printout`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Determines whether identified putative TF-binding sites for each species will be reported as a CSV file.


