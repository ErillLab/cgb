"""Module containing functions using NCBI Entrez utility.

See NCBI Entrez page (http://www.ncbi.nlm.nih.gov/books/NBK3837/) and Biopython
(http://biopython.org/DIST/docs/tutorial/Tutorial.html) tutorial for more
information.
"""

import os

from Bio import Entrez
from misc import directory
from my_logger import my_logger
import time

sleep_time = 0

# The directory used to save NCBI records for later use.
ENTREZ_DIRECTORY = directory('entrez_cache')


def set_entrez_email(email_address):
    Entrez.email = email_address

def set_entrez_apikey(api_key):
    Entrez.api_key = api_key

def set_entrez_delay(delay):
    global sleep_time
    sleep_time = delay

def get_genome_record(accession):
    """Gets the genome record from NCBI RefSeq."""
    global sleep_time
    genbank_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    if not os.path.isfile(genbank_file):
        # Download and save Genbank record
        my_logger.info("Downloading %s" % accession)
        handle = Entrez.efetch(db='nuccore', id=accession,
                               rettype='gbwithparts', retmode='text')
        record = handle.read()
        # add further delay if NCBI traffic is high and 
        # "HTTP Error 429: Too Many Requests" is received
        time.sleep(sleep_time)
        with open(genbank_file, 'w') as f:
            f.write(record)

    handle = open(genbank_file)
    return handle.read()


#takes an accession number for the protein, gets the record from NCBI
#(unless it is already stored locally) and saves it to file locally
#returns the object in the local file
def get_protein_record(accession):
    """Fetches the protein record from NCBI Protein database."""
    global sleep_time
    protein_file = os.path.join(ENTREZ_DIRECTORY, accession+'.gb')
    #if file not locally available, fetch and save locally (in ENTREZ_DIRECTORY cache)
    if not os.path.isfile(protein_file):
        # Download and save file
        handle = Entrez.efetch(db='protein', id=accession,
                               rettype='gb', retmode='text')
        record = handle.read()
        # add further delay if NCBI traffic is high and 
        # "HTTP Error 429: Too Many Requests" is received
        time.sleep(sleep_time)
        with open(protein_file, 'w') as f:
            f.write(record)

    #read file and return object
    handle = open(protein_file)
    return handle.read()
