from Bio import Entrez

Entrez.email = 'sefa1@umbc.edu'


def get_genome_record(accession):
    """Gets the genome record from NCBI RefSeq."""
    handle = Entrez.efetch(db='nuccore', id=accession,
                           retmode='gbwithparts', rettype='text')
    return handle.read()


def get_protein_record(accession):
    """Fetches the protein record from NCBI Protein database."""
    handle = Entrez.efetch(db='protein', id=accession,
                           rettype='gb', retmode='text')
    return handle.read()
