"""Class for performing BLAST."""

import os
import logging

from Bio.Blast import NCBIXML

from misc import temp_file_name


class BLAST:
    def __init__(self, seq_fasta, db_type, prefix=""):
        """Creates a BLAST database from the given fasta string.

        Args:
            seq_fasta (string): collection of sequences in FASTA format.
            db_type (string): 'nucl' or 'prot'
        """
        prefix = prefix + '_'
        self._in_file = temp_file_name(prefix=prefix, suffix='_input.fasta')
        self._db_file = temp_file_name(prefix=prefix, suffix='_blast.db')
        self._log_file = temp_file_name(prefix=prefix, suffix='_makeblastdb.log')
        self._seq_fasta = seq_fasta
        self._db_type = db_type
        self.makeblastdb()      # Create BLAST database.

    def makeblastdb(self):
        """Creates a BLAST database."""
        with open(self._in_file, 'w') as f:
            f.write(self._seq_fasta)
        cmd = 'makeblastdb -in {inp} -out {out} -logfile {log} -dbtype {db_type}'.format(
            inp=self._in_file, out=self._db_file, log=self._log_file,
            db_type=self._db_type)
        logging.debug(cmd)
        os.system(cmd)

    def search(self, blast_program, query, eval):
        """Runs BLAST to search query sequence in the target database.

        Args:
            blast_program (string): BLAST flavor to run. 'tblastx' or 'tblastn'
            query (string): query sequence in FASTA format
            eval (float): E-value threshold
        Returns:
            Bio.Blast.Record.Blast object
        """
        assert blast_program in ['tblastn', 'tblastx']
        output_file = temp_file_name()
        query_file = temp_file_name()
        with open(query_file, 'w') as f:
            f.write(query)
        cmd = '{prog} -query {q} -db {db} -evalue {e} -out {out} -outfmt 5'.format(
            prog=blast_program, q=query_file, db=self._db_file, e=eval,
            out=output_file)
        logging.debug(cmd)
        os.system(cmd)

        # Parse results
        with open(output_file) as results_handle:
            blast_record = NCBIXML.read(results_handle)
        return blast_record

    def tblastx(self, query, eval=10**-3):
        """Runs tblastx."""
        return self.search('tblastx', query, eval)

    def tblastn(self, query, eval=10**-3):
        """Runs tblastn."""
        return self.search('tblastn', query, eval)

    @staticmethod
    def get_best_hit(blast_record):
        """Returns the locus_tag of the best BLAST hit."""
        try:
            return blast_record.alignments[0].hit_def
        except IndexError:
            # no blast hits
            raise BlastNoHitFoundException

    @staticmethod
    def get_e_value(blast_record):
        """Returns the e-value of the best BLAST hit."""
        return blast_record.descriptions[0].e


class BlastNoHitFoundException(Exception):
    pass
