import os

from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'sefa1@umbc.edu'

ENTREZ_DIRECTORY = '/Users/sefa/Downloads/'


def get_genome_record(accession):
    """Gets the genome record from NCBI RefSeq."""
    genbank_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    if not os.path.isfile(genbank_file):
        # Download and save Genbank record
        handle = Entrez.efetch(db='nuccore', id=accession,
                               retmode='gbwithparts', rettype='text')
        with open(genbank_file, 'w') as f:
            f.write(handle.read())

    handle = open(genbank_file)
    return handle.read()



def get_protein_record(accession):
    """Fetches the protein record from NCBI Protein database."""
    handle = Entrez.efetch(db='protein', id=accession,
                           rettype='gb', retmode='text')
    return handle.read()
