===
cgb
===

``cgb`` is a Python library for comparative genomics of transcriptional
regulation in Bacteria.

.. sectnum::


About
-----

Given binding site evidence from one or more reference organisms, and a set of
target genomes of interest, ``cgb`` can be used to

- predict operons in each target genome
- identify putative binding sites on promoter regions of each target genome
- compute the posterior probability of regulation of each operon
- detect orthologs between organisms and create orthologous groups of genes
- perform ancestral state reconstruction on orthologous groups to analyze the
  evolution of transcriptional regulation of each gene


Get started
-----------

.. code-block:: python

   import cgb
   json_input_file = 'test_input.json'  # See below for the format
   cgb.go(json_input_file)

Input
-----

``cgb`` expects the input in JSON format. Below is a sample input file followed
by descriptions for each field.


.. code-block:: json

   {
    "TF": "LexA",
    "motifs": [
        {
            "protein_accession": "NP_217236.2",
            "sites": [
                "AAATCGAACATGTGTTCGAGTA",
                "GTCTCGAACATGTGTTCGAGAA",
                "GTATCGAACAATTGTTCGATAT",
                "GAATCAAACATGTGTTCGACAG",
                "TATTCGAACATGTATTCGAGTA"
            ]
        },
        {
            "protein_accession": "WP_003857389.1",
            "sites": [
                "TATGCGAACGTTTTTTCTAAAT",
                "TGATCGCAATTGTGTGCTAAAA",
                "TATTAAAACACTTGTTCTAAAC",
                "TAGTCGAACATGTGAACGGTAT",
                "AATACTGACAGAGGTTCGAATA",
                "ATCTCGAACACTCGTACCATTT",
                "ATTTCGAACAGTTGTGCGTGTA",
                "TATTCGAAAACTTTTCCGATCA",
                "TCCTCAAAAAAGTGGTCTAATG"
            ]
        }
    ],

    "genomes": [
        {
            "name": "ace",
            "accession_numbers": ["NC_008578.1"]
        },
        {
            "name": "cgl",
            "accession_numbers": ["NC_003450.3"]
        },
        {
            "name": "lxy",
            "accession_numbers": ["NC_006087.1"]
        }
    ],

    "prior_regulation_probability": 0.03,
    "phylogenetic_weighting": true,
    "site_count_weighting": true,
    "posterior_probability_threshold": 0.5
    }

Two mandatory input parameters are the list of reference motifs and target
genomes.

- The field ``motifs`` contains one or more motifs. Each motif is described by
  two sub-fields: ``protein_accession`` and ``sites``.

- The ``genomes`` field contains the list of target genomes to be used in the
  analysis. Each genome is described by two fields: ``name`` and
  ``accession_numbers``. The field ``accession_numbers`` could have multiple
  accession numbers, one for each chromosome/plasmid.

Other input parameters are optional.

- ``prior_regulation_probability``, the prior probability of regulation. Used
  by Bayesian estimation of probability of regulation.
- ``phylogenetic_weighting``. If true, the binding evidence from multiple
  reference organisms are weighted according to their phylogenetic distances to
  each target genome.
- ``site_count_weighting``. If true, the binding evidence from each reference
  organism is weighted by the binding site collection size.
- ``posterior_probability_threshold``. The genes/operons with posterior
  probability of regulation less than provided value are not reported.


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
