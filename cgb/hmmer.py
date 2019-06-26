"""Library for performing HMMER searches"""

import os
import subprocess

from math import log
from .my_logger import my_logger

from Bio import SearchIO

from Bio import SeqIO

HMMIN_FILENAME='/tmp/hmm_input.fasta'
HMMOUT_FILENAME='/tmp/hmm_tblout.tab'

def call_hmmscan(hmmscan_args):
    """Calls hmmersearch with given parameters

    Args:
        hmmscan_args (string): command line arguments, including query file
        and hmmer database file.
    """
    subprocess.call(hmmscan_args, stdout=open(os.devnull,'w'))

def process_hmmscan():
    """Processes the tab file resulting from a hmmscan search.
    """
    ###parse and return the result using hmmer3-tab format
    try:
        result=SearchIO.read(HMMOUT_FILENAME, 'hmmer3-tab')
        return(result)
    except ValueError:
        return([])

def run_hmmscan(prot_seq,eval,dbpathname):
    """Runs hmmscan search using provided protein sequence, the user
       configuration arguments for the search, and the specified mode.

       Args:
           prot_seq (seq object): sequence to be searched against HMM database
           eval (float): the e-value to use in hmmer
           dbpathname (string): the PATH (including filename) to HMMER database
    """

    #write to file the sequence in FASTA format
    with open("/tmp/hmm_input.fasta", "w") as file_handle:
        SeqIO.write(prot_seq, file_handle, "fasta")

    #compose the command line for hmmscan using the parameters in the user
    #input file: the hhmmer e-value limit and the hmmer database path and name
    hmmscan_args = ['hmmscan','--tblout', HMMOUT_FILENAME, '-E', eval, \
                    dbpathname, HMMIN_FILENAME]

    #call hmmscan with command line arguments
    call_hmmscan(hmmscan_args)

def run_COG_hmmscan(prot_seq,usr_inp):
    """Calls the run_hmmscan function with the appropriate database
    """
    #assign basic hmmscan parameters
    eval=str(usr_inp.hmmer_eval)

    #obtain a valid path for a HMMER database
    dbpathname = None
    try:
            dbpathname=usr_inp.COG_dbname

    except ValueError:
        my_logger.\
        info('***WARNING: Attempting COG HMMER query without defined database name')
        exit(0)

    if dbpathname == None:
        my_logger\
        .info('***WARNING: Attempting COG HMMER query without defined target database')
        exit(0)

    #call the hmmscan function with inferred evalue and db file path
    run_hmmscan(prot_seq,eval,dbpathname)
    
def run_eggNOG_hmmscan(prot_seq,usr_inp):
    """Calls the run_hmmscan function with the appropriate database
    """
    #assign basic hmmscan parameters
    eval=str(usr_inp.hmmer_eval)

    #obtain a valid path for a HMMER database
    dbpathname = None
    try:
            dbpathname=usr_inp.eggNOG_dbname

    except ValueError:
        my_logger.\
        info('***WARNING: Attempting NOG HMMER query without defined database name')
        exit(0)

    if dbpathname == None:
        my_logger\
        .info('***WARNING: Attempting NOG HMMER query without defined target database')
        exit(0)

    #call the hmmscan function with inferred evalue and db file path
    run_hmmscan(prot_seq,eval,dbpathname)

def run_PFAM_hmmscan(prot_seq,usr_inp):
    """Calls the run_hmmscan function with the appropriate database
    """
    #assign basic hmmscan parameters
    eval=str(usr_inp.hmmer_eval)

    #obtain a valid path for a HMMER database
    dbpathname = None
    try:
            dbpathname=usr_inp.PFAM_dbname

    except ValueError:
        my_logger.\
        info('***WARNING: Attempting PFAM HMMER query without defined database name')
        exit(0)

    if dbpathname == None:
        my_logger\
        .info('***WARNING: Attempting PFAM HMMER query without defined target database')
        exit(0)

    #call the hmmscan function with inferred evalue and db file path
    run_hmmscan(prot_seq,eval,dbpathname)

def process_COG_hmmscan(usr_inp):
    """Calls the process_hmmscan function to analyze hmmscan output and
       postprocesses the result according to the options in cofiguration file.
       For COG queries, the user can specify a maximum number of results,
       and/or a limiting exponent jump in e-value. For instance, if limiting
       e-value jump is set to 3, and the first COG has e-value 10e-30, the
       function will only report back up to 10e-27.
    """

    #call the hmmscan processing function, which returns a list of hit objects
    COGs=process_hmmscan()

    COG_list=[]
    mineval=0.0

    #determine if there is at least a hit
    if len(COGs)>0:
        #if there is a hit, get its e-value
        mineval=COGs[0].evalue
        logmineval=log(mineval+1e-250,10)

        cnt=1
        #for every result
        for res in COGs:
            #check if e-value is below jump
            if ((log(res.evalue+1e-250,10) - logmineval) \
                 < usr_inp.OGejump):
                #check if we are still below max number of COGs to report
                if cnt<=usr_inp.maxNOG:
                    element={'ID' : res.id, 'eval' : res.evalue,\
                             'desc' : res.description}
                    #append COG ID
                    COG_list.append(element)
            cnt=cnt+1

    #return COG_list
    return COG_list

def process_eggNOG_hmmscan(usr_inp):
    """Calls the process_hmmscan function to analyze hmmscan output and
       postprocesses the result according to the options in cofiguration file.
       For eggNOG queries, the user can specify a maximum number of results,
       and/or a limiting exponent jump in e-value. For instance, if limiting
       e-value jump is set to 3, and the first NOG has e-value 10e-30, the
       function will only report back up to 10e-27.
    """

    #call the hmmscan processing function, which returns a list of hit objects
    NOGs=process_hmmscan()

    NOG_list=[]
    mineval=0.0

    #determine if there is at least a hit
    if len(NOGs)>0:
        #if there is a hit, get its e-value
        mineval=NOGs[0].evalue
        logmineval=log(mineval+1e-250,10)

        cnt=1
        #for every result
        for res in NOGs:
            #check if e-value is below jump
            if ((log(res.evalue+1e-250,10) - logmineval) \
                 < usr_inp.OGejump):
                #check if we are still below max number of NOGs to report
                if cnt<=usr_inp.maxNOG:
                    #clean up and append NOG ID
                    #assumes that NOG ID is of the form:
                    #'bctoNOG.ENOG4109EA4.meta_raw'
                    #and we want to keep 'ENOG4109EA4'
                    NOG=res.id
                    st = NOG.find('ENOG')
                    ed = NOG.find('.',st)
                    element={'ID' : NOG[st:ed],  'eval' : res.evalue,\
                             'desc' : res.description}                    
                    NOG_list.append(element)
            cnt=cnt+1

    #return NOG_list
    return NOG_list


def process_PFAM_hmmscan(usr_inp):
    """Calls the process_hmmscan function to analyze hmmscan output and
       postprocesses the result according to the options in cofiguration file.
       For PFAM queries, the user can specify a maximum number of results.
    """

    #call the hmmscan processing function, which returns a list of hit objects
    PFAMs=process_hmmscan()

    PFAM_list=[]

    #determine if there is at least a hit
    if len(PFAMs)>0:
        cnt=1
        #for every result
        for res in PFAMs:
            #check if we are still below max number of PFAMs to report
            if cnt<=usr_inp.maxPFAM:
                #clean up assign PFAM ID accession
                #returned PFAM accession has version number, which is dropped
                #e.g. PFAM3213.3 --> PFAM3213
                PFAM=res.accession
                element={'ID' : PFAM.split('.')[0], 'eval' : res.evalue,\
                         'desc' : '//'.join([res.id,res.description])}                
                PFAM_list.append(element)
            cnt=cnt+1

    #return PFAM_list
    return PFAM_list