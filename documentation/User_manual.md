CGB user manual
===============

Introduction
------------

### What is CGB?

CGB is a client-based application for comparative genomics analyses of transcriptional regulatory networks in Bacteria (and Archaea). Using multiple sources of information on the reported binding specificity of a transcription
factor (TF), CGB scans bacterial genomes looking for putative regulated operons, assesses their probability of regulation and analyzes its conservation across species using ancestral state reconstruction, providing an easily interpretable overview of a transcriptional regulatory network across multiple bacteria.

### Why CGB?

Genomics data is probably the most abundant resource in microbial bioinformatics, and it can be effectively leveraged to reconstruct transcriptional regulatory networks across multiple bacterial genomes, providing an evolutionary outline of the conservation and changes in a bacterial transcriptional regulatory network. A few tools for comparative genomics of regulatory networks in Bacteria are available, but they are server-oriented, limiting the extent of the analysis to preprocessed genomes and constraining the fine-tuning of analysis parameters and the availability of raw output data for analysis.

CGB provides a fully-customizable, open-source scalable solution for comparative genomics of transcriptional regulatory networks in Bacteria. By applying well-justified heuristics, CGB speeds up computation time, enabling customized analyses to be performed on the client side. Furthermore, CGB uses a rigorous, generalizable Bayesian framework to compute regulation probabilities and reconstruct ancestral regulation state, enabling users to draw informed conclusions on the evolutionary history of a regulatory network and the likelihood and relevance of the regulation of certain elements therein.

### Main features

CGB introduces a plethora of innovative features that enable fast, personalized and easily interpretable comparative analyses of transcriptional regulatory networks in Bacteria.

-   Open source code, fully customizable

-   Client based, runs locally on a modest desktop setup

-   Can use any genome (complete or incomplete) available on the NCBI repository

-   Automatic handling of composite (multiple chromosomes, plasmids) genomes

-   Allows integrating multiple sources of evidence on TF-binding specificity

-   Generalizable TF-binding model

-   Integrated operon prediction with smart operon-splitting

-   Computes easily interpretable Bayesian posterior probability of regulation

-   Gene-based integration of analysis results

-   Integrated ancestral state reconstruction for gene regulation

-   Graphical and complete raw data (csv) output on all aspects of the analysis

Basic operation
---------------

### Dependencies

Apart from the python 2.7 libraries listed on ```requirements.txt```, CGB makes use of three external programs:

-   NCBI BLAST

-   ClustalO

-   BayesTraits

### Outline

CGB takes as input three main elements: (1) a set of genomes to be analyzed and (2) a set of proteins (transcription factor instances) for which (3) experimental evidence of binding (a binding site collection) is available. CGB then integrates the available experimental data into TF-binding models adapted to each target species under analysis, and uses these TF-binding models to analyze the probability of regulation of each gene. It then predicts operons in each genome and splits them if an internal gene is likely regulated. Orthologs for likely-regulated genes are detected across all target genomes, and the evidence of regulation for any given gene is assessed using ancestral state reconstruction on the posterior probability of regulation of the operon harboring the gene in each species.

### Pipeline

CGB implements a fairly linear comparative genomics pipeline that can be summarized as follows:

-   Input processing

-   Protein & genome fetch and instantiation

-   Automated detection of TF in target genomes

-   TF-based phylogeny generation

-   Site collection instantiation (experimental data)

-   Computation of evidence transfer weights

-   Generation of species-specific TF-binding models

-   Prior estimation

-   Gene-wise analysis of posterior probability of regulation

-   Operon prediction and split

-	Regulon inference and ortholog prediction

-   Ancestral state reconstruction

-   Output generation

Pipeline rationale and details
------------------------------

### Assumptions

CGB assumes that on the species under analysis:

-   The TF is conserved

-   The TF-binding specificity is conserved to some extent [i.e. shared compatible motif]

-   There are no close paralogs of the TF or, if there are, they share a conserved motif

CGB will run if these assumptions do not hold, but caution should be applied when analyzing genomes in which the user suspects that TF gene duplication has occurred.

### Input processing

To facilitate interoperability, all required inputs are encoded in a JSON file (as lists of dictionaries). The JSON file input is saved onto a UserInput Python object, which essentially provides get methods.

### Protein & genome fetch and instantiation

After reading the JSON input, the pipeline instantiates the protein and genome records. For all TF instances associated with experimental support, the protein record for the accession is fetched from NCBI. For all the chromids (chromosomes or plasmids) listed as components of a target genome, their NCBI GenBank record is retrieved and integrated into a ```genome``` object.

### Automated detection of TF in target genomes

Using the provided TF instances with experimental support, a BLAST search with 10\^{-3} e-value limit is performed against each target genome. The best first hit (from any reference TF) is considered the TF homolog in each genome [Note: this assumes no close TF paralogs are present or relevant]. Genomes not containing a TF homolog are dropped from the analysis.

### TF-based phylogeny generation

A phylogeny of the species under analysis (and species with reference information) is generated using the protein sequence of the TF homologs detected in each target genome. This phylogeny is used to perform ancestral state reconstruction and to weight the contribution of each reference binding site collection to create species-specific TF-binding models and set species-specific priors for regulation.

The rationale for choosing the TF as the element to perform phylogenetic inference on is two-fold. On the one hand, it is by construction the only gene known to be preserved in all species under analysis. On the other hand, for the purposes of inferring TF-binding models and priors, and making inferences on regulatory network evolution, it is the most relevant protein. While other phylogenetic makers can conceivably be used, their presence in all genomes is not guaranteed and the relevance of their particular phylogenetic inference to the evolution of the regulatory network is not necessarily warranted.

### Site collection instantiation (experimental data)

Every reference TF provided in the input file is associated with a collection of *aligned* experimental TF-binding sites. For each TF, the site list is parsed and associated with the TF in a ```collection``` object.

### Computation of evidence transfer weights

The first step in the main pipeline loop is to compute the weighting scheme to transfer information into the genome. This information is, essentially, the TF-binding model and the prior for regulation used in Bayesian estimation. Weights are multiplicative factors on the additive normalized sum of the contributions to the TF-model or prior probability for that genome.

Experimental binding site information from reference genomes can be transferred to a given target genome in three different ways:

-   Unweighted: all information is transferred using a simple arithmetic mean (all weights set to 1)

-   Phylogenetic weighting: information is transferred using a weighted mean, with phylogenetic distance *decreasing* the contribution of a reference to a target

-   Phylogenetic weighting with site-count weights: information is transferred using a weighted mean that combines (as a product) the phylogenetic and site count weights; the number of sites in a given collection *increases* the contribution of that reference to target genomes, modulated by the phylogenetic distance

The estimation of priors for regulation can use the unweighted approach or phylogenetic weighting.

#### Phylogenetic weighting

Phylogenetic weighting is performed using a "*clustalesque*" weighting scheme. The idea here is the same as in CLUSTALW. If a branch of the tree has a lot of experimental evidence, and another only a bit, the "populated" branch is going to dominate, and set all motifs and priors according to that side of the tree (even in spite of contradictory evidence somewhere else on the tree). To counteract this, the *clustalesque* approach dilutes the contribution of a populated branch by recomputing the distance of the elements in that branch to the target genome, essentially making them "more distant" (proportionally to how many informative elements lie below the branch).

If *clustalesque* weighting is activated, the TF-based tree is first re-rooted using the genome of interest as the outgroup. The branch length of each clade (what connects it to the upper part of the tree) is multiplied by the number of informative elements within the clade (this pushes the clade away proportionally to how much information-bearing species it contains) and this is done for all clades.

Weights are then computed as follows. The target genome is identified on the tree, and for each site collection, we find its corresponding node on the tree and measure the distance (*clustalesque* modified or not) between the target genome and the site collection node. These distances are then transformed to weights, by computing [1.0 - (dist - min(distances))/max(distances). After all weights have been computed for that target species, weights are normalized to add up to one.

The use of (1-(dist-min))/max to estimate weights from distances seeks to get a 1 (maximal) contribution if the distance between target and reference is small (min), or a very small contribution (max-min)/max if the reference is very far away (max).

#### Site count weighting

For site count weighting, the weights are literally the number of sites available in each reference collection, normalized (after all reference collections have been assessed) to add up to 1. This allows larger collections to dominate (i.e contribute more to targets), but still takes into account the contribution of small collections if they are close to the target.

The rationale behind this approach is that we are considering site collections as "approximations" to TF-binding models. Ideally, site collections would be complete. In such a case, if a TF binds only one site in the genome, that is essentially your binding model. However, in real life most collections are going to be incomplete (especially those with only one site), and you don't want the incomplete model to dominate its phylogenetic neighborhood. With the non-absurd assumption that the number of genes regulated by a TF will be roughly similar (same order-of-magnitude) across species, using the number of sites in the reference collections to weight up or down their contribution to targets (modulated by phylogeny) corrects naturally for the assumption that provided models are complete (i.e. it assumes that smaller collection are most likely due to experimental undersampling).

Phylogenetic weighting and site count weighting weights are combined by multiplying them (for each reference, multiply its phylogenetic and site-count weight to target).

### Generation of species-specific TF-binding models

Once weights have been computed (or with the default weights set to one), the pipeline computes the TF-binding model (PSSM model) of each target genome. These models are computed by obtaining the PSWM (Position-Specific Weight Matrix) stemming from each reference collection, multiplying the frequencies for each base in each column with the computed weights between the reference and the target, and creating the target PSWM as the weighted average of each frequency. Once the target PSWM is generated, the PSSM model is computed (using homogeneous 0.25 base frequencies by default) and its background and foreground score distribution is estimated on a randomly selected set of promoter regions from the target genome and the PSWM itself, respectively. These two score distributions are used in the computation of the posterior probability of binding to a given sequence.

#### Statistical model
Given a TF-binding model (e.g. a PSWM) we can generate a scoring system (e.g. PSSM) that defines the likelihood of the TF binding a given sequence. To take into account that the TF may bind a given DNA sequence on either strand (in fact, making contacts with both strands), we define:

```PSSM(s_i)=log_2\left ( 2^{PSSM()s_i^f)} + 2^{PSSM(s_i^r)} \right )```
<a href="https://www.codecogs.com/eqnedit.php?latex=PSSM(s_i)=log_2\left&space;(&space;2^{PSSM()s_i^f)}&space;&plus;&space;2^{PSSM(s_i^r)}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?PSSM(s_i)=log_2\left&space;(&space;2^{PSSM()s_i^f)}&space;&plus;&space;2^{PSSM(s_i^r)}&space;\right&space;)" title="PSSM(s_i)=log_2\left ( 2^{PSSM()s_i^f)} + 2^{PSSM(s_i^r)} \right )" /></a>

Given this scoring function, we can define the following:

*Expected distribution of scores on non-regulated promoters*
```B\propto N(\mu_g,\sigma_g^2)```
<a href="https://www.codecogs.com/eqnedit.php?latex=B\propto&space;N(\mu_g,\sigma_g^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B\propto&space;N(\mu_g,\sigma_g^2)" title="B\propto N(\mu_g,\sigma_g^2)" /></a>

*Expected distribution of scores on regulated promoters*
```R\propto \alpha N(\mu_m,\sigma_m^2) + (1-\alpha)N(\mu_g,\sigma_g^2)```
<a href="https://www.codecogs.com/eqnedit.php?latex=R\propto&space;\alpha&space;N(\mu_m,\sigma_m^2)&space;&plus;&space;(1-\alpha)N(\mu_g,\sigma_g^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R\propto&space;\alpha&space;N(\mu_m,\sigma_m^2)&space;&plus;&space;(1-\alpha)N(\mu_g,\sigma_g^2)" title="R\propto \alpha N(\mu_m,\sigma_m^2) + (1-\alpha)N(\mu_g,\sigma_g^2)" /></a>

where ```\alpha``` defines a mixing coefficient between the distribution expected among "true" binding sequences (*m*) and the distribution among non-sites (*g*). Both distributions are approximated by a normal, parametrized by the data (the PSWM and the genome, respectively).

### Prior estimation

The user can define an analysis-wide prior for regulation. The prior for regulation is used on the computation of the probability of regulation of any given sequence. Essentially, we seek to compute P(R\|S)=P(S\|R)·P(R)/(P(S\|R)·P(R)+P(S\|B)·P(B)).

``P(R|D)=\frac{P(D|R)P(R)}{P(D)}=\frac{P(D|R)P(R)}{P(D|R)P(R)+P(D|B)P(B)}`` `
<a href="https://www.codecogs.com/eqnedit.php?latex=P(R|D)=\frac{P(D|R)P(R)}{P(D)}=\frac{P(D|R)P(R)}{P(D|R)P(R)&plus;P(D|B)P(B)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P(R|D)=\frac{P(D|R)P(R)}{P(D)}=\frac{P(D|R)P(R)}{P(D|R)P(R)&plus;P(D|B)P(B)}" title="P(R|D)=\frac{P(D|R)P(R)}{P(D)}=\frac{P(D|R)P(R)}{P(D|R)P(R)+P(D|B)P(B)}" /></a>

If we assume positional independency, the likelihoods can be estimated from the background/foreground score distributions as follows:
```P(D|R)=\prod_{s_i\in D} L(s_i|R)= \prod_{s_i\in D} L\left (s_i | \alpha N(\mu_m,\sigma_m^2) + (1-\alpha)N(\mu_g,\sigma_g^2)  \right )```
<a href="https://www.codecogs.com/eqnedit.php?latex=P(D|R)=\prod_{s_i\in&space;D}&space;L(s_i|R)=&space;\prod_{s_i\in&space;D}&space;L\left&space;(s_i&space;|&space;\alpha&space;N(\mu_m,\sigma_m^2)&space;&plus;&space;(1-\alpha)N(\mu_g,\sigma_g^2)&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P(D|R)=\prod_{s_i\in&space;D}&space;L(s_i|R)=&space;\prod_{s_i\in&space;D}&space;L\left&space;(s_i&space;|&space;\alpha&space;N(\mu_m,\sigma_m^2)&space;&plus;&space;(1-\alpha)N(\mu_g,\sigma_g^2)&space;\right&space;)" title="P(D|R)=\prod_{s_i\in D} L(s_i|R)= \prod_{s_i\in D} L\left (s_i | \alpha N(\mu_m,\sigma_m^2) + (1-\alpha)N(\mu_g,\sigma_g^2) \right )" /></a>

The priors P(R) and P(B) represent our belief that an operon might (or not) be regulated a priori. P(B) is by definition 1-P(R). If an overarching P(R) prior is not provided by the user, the pipeline estimates the prior for each reference species as the ratio between reported binding sites (assuming one site per regulatory region) and the number of operons in a species. Thus, if 10 sites are reported for a TF in a species with a 1,000 operons, the prior for regulation on any given operon is P(R)=10/1,000=0.01, and P(B)=0.99.

The priors for all different reference species are first computed and then merged into priors for target genomes through the same weighting approach used to generate the TF-binding models, but including only the phylogenetic component (not site counts; XXX verify). It is worth noting that this method will tend to underestimate the prior probability of regulation, because reported collections are not complete. This will in turn bias the inferred posterior probability of regulation P(R\|S), making it less responsive (XXX check?). The user can adjust the inferred ratios globally through a tuning parameter (XXX; ToDo)

### Gene-wise analysis of posterior probability of regulation

Once TF-binding models have been instantiated in each target genome, the pipeline performs an analysis of posterior probability of regulation for each *gene*, computed over its upstream region. This posterior is stored in the ```gene``` object.

### Operon prediction and split

Operon prediction is performed using an intraoperon intergenic distance threshold. That is, two neighboring genes in the same strand are considered to be in the same operon if their intergenic distance is less than or equal to the intraoperon intergenic distance threshold.

Because different genomes show different degrees of compaction, CGB estimates the intraoperon intergenic distance threshold on each genome. It does so as follows. First, all *directons* in the genome are identified. Directons are consecutive genes in the same strand. The pipeline then computes the mean intraoperon intergenic distance between the first two genes in any directon that has contains than one gene and sits opposed to another directon [(2<-(x)<-1<- ->1->(x)->2)]. These directons are choosen to estimate the mean intraoperon intergenic distance because, without further knowledge, they are the most likely bona fide operons (i.e. they are in opposing strands and hence the first operon gene must be the start of a genuine operon; the assumption then is that for a substantial number of cases the operon will contain at least two genes and hence the mean intergenic distance between first and second gene is a decent estimate of intraoperon intergenic distance).

The mean intraoperon intergenic distance is uses as the intraoperon intergenic distance threshold.

#### Operon split

After operons has been predicted, the pipeline assesses whether any operons should be split. 
Given that functional binding sites are very unlikely to occur anywhere except in regulated promoters, the pipeline analyzes whether the upstream region of any gene contains a high-scoring putative site. Pperons that contain a high-scoring putative site upstream of an internal gene (i.e. any gene of the operon except the first one) are split into two operons at the point where the putative binding site lies. It should be noted that the PSSM score threshold that determines what is considered a high-scoring TF-binding site is used only to improve operon prediction, not to infer regulation. The selected PSSM score threshold satisfies the equality between the negative logarithm of the false positive rate (FPR) and the information content (IC) of the motif (Hertz and Stormo, 1999).


This intraoperon intergenic distance threshold can be scaled up or down by a ```\sigma``` user parameter. The rationale for this adjustment is that operon splitting lowers the chances of *missing* regulated genes by packing them within an "incorrect" operon structure (that shows no evidence of regulation). When predicting operons for regulatory analysis, there are two confounding problems: (1) incorrectly globbling up genes (making operons larger than they ought to be) and (2) incorrectly splitting operons. The use of an operon splitting strategy based on the presence of likely TF-binding sites upstream of *internal* genes helps to palliate the gobbling problem (1). This allows using larger intergenic distance thresholds (more gobbling-prone), which decreases the chances of inadvertently splitting a bona fide operon (2).

### Regulon inference and ortholog prediction

After predicting and splitting operons, all operons with a posterior probability of regulation larger than the user-specified threshold (0.5 by default) are considered, in each genome, members of the putative regulon. To enhance runtime, orthologs are only computed for members of the putative regulon in each genome. Orthologs are computed genome-wise using the best-reciprocal BLAST hit principle. Orthologs are stored in the ```OrthologousGroup``` object.

### Ancestral state reconstruction

Once orthologs have been detected, the ancestral state of the posterior probability of regulation for each gene is reconstructed using BayesTraits, by mapping for each ortholog the posterior probability of regulation on the phylogenetic tree generated with the TF protein sequence. Since BayesTraits performs ancestral reconstruction on discrete characters, CGB generates bootstrap replicates in which the posterior is discretized into a binary state (regulated or not), based on the posterior value (i.e. a 0.7 posterior results in 70% of replicates assigning that gene as regulated). A third state is used to designate the case in which the species does not contain an ortholog for the gene (which is not the same as inferring that it is not regulated). The reconstructed states for all bootstrap replicates are then mapped back into probabilities using the frequency of each reconstructed state as the MLE for the posterior.

### Output generation

Input
-----

-   motifs: a number of experimentally-determined TF-binding sites, together
    with the protein accession for the TF and the species name (used for the
    tree visualization)

-   genomes: all the genomes to be analyzed, composed of a name and a list of
    accession numbers (can be draft/incomplete genomes)

-   parameters:

    -   prior\_regulation\_probability: genomic prior for *operon* regulation

    -   phylogenetic\_weighting: ON/OFF

    -   site\_count\_weighting: ON/OFF

    -   posterior\_probability\_threshold: posterior cut for regulation
        assignment

    -   posterior\_probability\_threshold\_for\_operon\_prediction: posterior
        cut for splitting operons

    -   ancestral\_state\_reconstruction: ON/OFF

    -   operon\_prediction\_distance\_tuning\_parameter: used to
        increase/decrease the computed mean intergenic distance

a BLAST search with 10\^{-3} TF homolog detection

modulator for prior estimation?

Number of sites from promoter regions used to estimate background?

Background frequencies for TF-model log-likelihood?

Define upstream region limits?

Output
------
``cgb`` saves all the output in the folder ``output`` created on the working
directory.

- ``user_PSWM/`` contains the user-provided binding motifs in JASPAR format.

- ``derived_PSWM/`` contains binding motifs in JASPAR format, tailored for each
  target genome combining all the evidence from each reference motif.

- ``identified_sites/`` contains identified binding sites and information such
  as their genomic locations, downstram regulated genes and their
  functions. Predicted binding site data is saved into CSV files, one for each
  target genome.

- ``operons/`` contains the operon predictions of each target genome, saved as
  CSV files.

- ``orthologs.csv`` contains the groups of orthologous genes and their
  probabilities of regulation.

- ``phylogeny.png`` is plot of the phylogenetic tree.

- ``ancestral_states.csv`` has the reconstructed state of each gene in all
  ancestral clades. For each target species and ancestral clades, the states
  are

  - ``P(1)``, the probability of TF *binding*
  - ``P(0)``, the probability of TF *not binding*
  - ``P(A)``, the probability of *absence* of the gene.

- ``plots/`` folder contains the visualization of the results.