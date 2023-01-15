#!/usr/bin/env python

"""
This module defines class BLAST for labas.

Dependencies: BLAST+, BioPython, Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 14 Jan 2023; the latest update: 15 Jan 2023.
"""

import os
import sys
import subprocess
from distutils import spawn
from collections import namedtuple

BLAST_parameters = namedtuple('BLAST_parameters', ['min_identity', 'min_qcov', 'max_evalue', 'max_hits'])

class BLAST:
    """ Execute BLAST tasks with the same query FASTA file """
    def __init__(self, query_fasta, min_identity, min_qcov, max_evalue, max_hits):
        """ Initiate private properties """
        self.__query_fasta = query_fasta  # Fasta file of query sequences
        self.__params = BLAST_parameters(str(min_identity), str(min_qcov), max_evalue, str(max_hits))
        return
    
    @property
    def query_fasta(self):
        return self.__query_fasta

    @property
    def params(self):
        return self.__params
    
    def search(self, subject_name, subject_fasta, outdir):
        """ Search query sequences against a database of the subject sequence using megaBLAST """
        output_tsv = os.path.join(outdir, subject_name + '__megaBLAST.tsv')
        header = ' '.join(['6', 'qseqid', 'sseqid', 'slen', 'pident', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',\
                           'sstart', 'send', 'sstrand', 'evalue', 'bitscore', 'sseq'])
        blastn_command = ['blastn', '-task', 'megablast', '-query', self.__query_fasta, '-subject', subject_fasta,\
                          '-out', output_tsv, '-perc_identity', self.__params.min_identity, '-qcov_hsp_perc', self.__params.min_qcov,\
                          '-evalue', self.__params.max_evalue, '-max_target_seqs', self.__params.max_hits, '-outfmt', header]
        process = subprocess.Popen(blastn_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        _, err = process.communicate()  # obtain outcomes (out and err) as two long strings (in bytes class though)
        #out = out.decode()
        err = err.decode()
        if len(err) > 0:
            print("\nblastn encountered an error:", file = sys.stderr)
            print(err, file = sys.stderr)
            sys.exit(1)
        return

    def check_executives(self):
        #makeblastdb_path = spawn.find_executable("makeblastdb")
        blastn_path = spawn.find_executable("blastn")
        #blast_installed = (makeblastdb_path != None and blastn_path != None)
        #if not blast_installed:
        if blastn_path == None:
            print("Error: could not find executives of BLAST.", file = sys.stderr)
            sys.exit(1)
        return

    """
    def make_blast_db(self, db_dir, subject_name, subject_fasta):
        ''' Create a BLAST nucleotide database in directory db_dir (not needed for a search between two FASTA files) '''
        db = os.path.join(db_dir, subject_name)
        makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', subject_fasta, '-out', db]
        process = subprocess.Popen(makeblastdb_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        _, err = process.communicate()  # Return bytes not strings
        err = err.decode()  # Convert bytes into a string
        if len(err) > 0:
            print("\nmakeblastdb encountered an error:", file = sys.stderr)
            print(err, file = sys.stderr)
            sys.exit(1)
        return db

    def clean(db):
        ''' Delete database files *.ndb, *.nhr, *.nin, etc. '''
        subprocess.run(f'rm {db}.*', shell = True)
        return
    """
