#!/usr/bin/env python

"""
This module defines class BLAST.

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
from lib.Hit import HIT_TABLE_COLS

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
        """
        Search query sequences against a database of the subject sequence using megaBLAST
        Return: content of a tsv file for BLAST results of the current sample
        """
        header = HIT_TABLE_COLS + ['sseq']
        blastn_command = ['blastn', '-task', 'megablast', '-query', self.__query_fasta, '-subject', subject_fasta,\
                          '-perc_identity', self.__params.min_identity, '-qcov_hsp_perc', self.__params.min_qcov,\
                          '-evalue', self.__params.max_evalue, '-max_target_seqs', self.__params.max_hits, '-outfmt', ' '.join(['6'] + header)]
        process = subprocess.Popen(blastn_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = process.communicate()  # Obtain outcomes (out and err) as two long strings (in bytes class though)
        out = out.decode()
        err = err.decode()
        if len(err) > 0:
            print("\nblastn encountered an error:", file = sys.stderr)
            print(err, file = sys.stderr)
            sys.exit(1)
        output_tsv = open(os.path.join(outdir, subject_name + '__megaBLAST.tsv'), 'w')  # Save raw BLAST outputs to subdirectory 1_blast
        if len(out) > 0:
            output_tsv.write('\t'.join(header) + "\n")
            output_tsv.write(out)
            out = out.splitlines()  # Convert lines into a list
        else:
            print(f"Warning (blast.search): no match of any query sequences was found in sample {subject_name} and therefore, no result file will be produced.", file = sys.stderr)
            out = None
        output_tsv.close()  # Leaves an empty file if len(out) = 0
        return out

    def check_executives(self):
        blastn_path = spawn.find_executable("blastn")
        if blastn_path == None:
            print("Error: could not find executives of BLAST.", file = sys.stderr)
            sys.exit(1)
        return
