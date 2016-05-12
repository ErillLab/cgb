cgb
===

``cgb`` is a Python library for comparative genomics of transcriptional
regulation in Bacteria.

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

Input format
------------

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

-
