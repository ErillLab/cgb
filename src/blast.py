"""Class for performing BLAST."""

import os

from Bio.Blast import NCBIXML


class BLAST:
    def __init__(self, seq_fasta, db_type):
        """Creates a BLAST database from the given fasta string.

        Args:
            seq_fasta (string): collection of sequences in FASTA format.
            db_type (string): 'nucl' or 'prot'
        """
        self._in_file = '/tmp/input.fasta'
        self._db_file = '/tmp/blast.db'
        self._log_file = '/tmp/makeblastdb.log'
        self._seq_fasta = seq_fasta
        self._db_type = db_type

    def makeblastdb(self):
        """Creates a BLAST database."""
        with open(self._in_file, 'w') as f:
            f.write(self._seq_fasta)
        cmd = 'makeblastdb -in {inp} -out {out} -logfile {log} -dbtype {db_type}'.format(
            inp=self._in_file, out=self._db_file, log=self._log_file,
            db_type=self._db_type)
        print cmd
        os.system(cmd)

    def tblastx(self, query, eval=10**-3):
        output_file = '/tmp/tblastx.out'
        query_file = '/tmp/query.fasta'
        with open(query_file, 'w') as f:
            f.write(query)
        cmd = 'tblastx -query {q} -db {db} -evalue {e} -out {out} -outfmt 5'.format(
            q=query_file, db=self._db_file, e=eval, out=output_file)
        print cmd
        os.system(cmd)

        # Parse results
        with open(output_file) as results_handle:
            blast_record = NCBIXML.read(results_handle)
        return blast_record
