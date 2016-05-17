"""Module containing functions using NCBI Entrez utility.

See NCBI Entrez page (http://www.ncbi.nlm.nih.gov/books/NBK3837/) and Biopython
(http://biopython.org/DIST/docs/tutorial/Tutorial.html) tutorial for more
information.
"""

import os

from Bio import Entrez
from misc import directory
from my_logger import my_logger

Entrez.email = 'sefa1@umbc.edu'

# The directory used to save NCBI records for later use.
ENTREZ_DIRECTORY = directory('entrez_cache')


def get_genome_record(accession):
    """Gets the genome record from NCBI RefSeq."""
    genbank_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    if not os.path.isfile(genbank_file):
        # Download and save Genbank record
        my_logger.info("Downloading %s" % accession)
        handle = Entrez.efetch(db='nuccore', id=accession,
                               retmode='gbwithparts', rettype='text')
        record = handle.read()
        with open(genbank_file, 'w') as f:
            f.write(record)

    handle = open(genbank_file)
    return handle.read()


def get_protein_record(accession):
    """Fetches the protein record from NCBI Protein database."""
    protein_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    if not os.path.isfile(protein_file):
        # Download and save file
        handle = Entrez.efetch(db='protein', id=accession,
                               rettype='gb', retmode='text')
        record = handle.read()
        with open(protein_file, 'w') as f:
            f.write(record)

    handle = open(protein_file)
    return handle.read()
