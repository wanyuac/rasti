#!/usr/bin/env python

"""
This module defines class Queries for labas.

Dependencies: Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 15 Jan 2023.
"""

from Bio import SeqIO

class Queries:
    """ Parse the FASTA file of query sequences and manage its sequence data """
    def __init__(self, fasta):
        self.__qlens = dict()
        qs = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        for q, s in qs.items():
            self.__qlens[q] = len(str(s.seq))
        return
    
    @property
    def query_names(self):
        return list(self.__qlens.keys())

    @property
    def query_num(self):
        return len(self.__qlens)

    def write_query_lengths(self, tsv):
        f = open(tsv, 'w')
        f.write('\t'.join(["Query", "Length_bp"]) + '\n')
        for q, l in self.__qlens.items():
            f.write('\t'.join([q, str(l)]) + '\n')
        f.close()
        return
