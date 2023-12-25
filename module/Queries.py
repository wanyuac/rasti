#!/usr/bin/env python

"""
This module defines class Queries.

Dependencies: Python 3

Copyright (C) 2023-2024 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 26 Dec 2023.
"""

from Bio import SeqIO
from collections import namedtuple

Query = namedtuple('Query', ['len', 'type'])  # Length (bp) and type (CDS, IS, etc)

class Queries:
    """
    This class handles the input multi-FASTA file of query sequences. It parses the FASTA file and
    manages its sequence data. This class does not store query sequences but their names and lengths.
    """
    def __init__(self, fasta):
        self.__queries = dict()  # sequence name : Query(sequence length, feature type)
        self.__cds = list()  # Names of coding sequences
        qs = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
        for q, s in qs.items():  # Element 's' is a SeqRecord.
            t = self.__extract_seq_type(s.description)  # Headers for CDSs: >seq CDS|annotations
            self.__queries[q] = Query(len = len(str(s.seq)), type = t)
            if t == 'CDS':
                self.__cds.append(q)
        return
    
    @property
    def query_names(self):
        return list(self.__queries.keys())

    @property
    def query_num(self):
        return len(self.__queries)
    
    @property
    def cds(self):
        """ Returns a list of the names of all query CDSs """
        return self.__cds
    
    @property
    def cds_num(self):
        """ Returns the number of CDSs """
        return len(self.__cds)

    def __extract_seq_type(self, h):
        descr = h.split(' ')[1]
        return descr.split('|')[0]

    def query_len(self, q):
        return self.__queries[q].len  # Value type: integer
    
    def query_type(self, q):
        return self.__queries[q].type  # Value type: string

    def write_query_lengths(self, tsv):
        """ Prints a summary TSV file of query sequences in Stage 2 of function detect """
        f = open(tsv, 'w')
        f.write('\t'.join(["Query", "Type", "Length"]) + '\n')
        for i, q in self.__queries.items():
            f.write('\t'.join([i, q.type, str(q.len)]) + '\n')
        f.close()
        return
