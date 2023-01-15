#!/usr/bin/env python

"""
This module defines class Hits for labas.

Dependencies: Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 15 Jan 2023.
"""

import os
import sys
from lib.Hit import Hit
from lib.Hit import HIT_TABLE_COLS

class Hit_tables:
    """ Manage raw BLAST outputs """
    def __init__(self):
        self.__hit_tables = dict()  # A dictionary {genome : table} of raw BLAST outputs
        self.__hit_counts = dict()  # Number of hits per sample
        return
    
    @property
    def table_num(self):  # Number of hit tables
        return len(self.__hit_tables)

    @property
    def sample_names(self):
        return list(self.__hit_tables.keys())

    def add_table(self, sample, hit_table):  # Parameter hit_table is a list of rows in BLAST's raw tabulated output.
        if hit_table != None:  # Each element of the dictionary __hit_tables stores information of rows in the raw output TSV file when the output is not empty. Otherwise, the table is None.
            self.__hit_counts[sample] = len(hit_table)  # Count the number of rows in the current output table of BLAST.
            ht = dict()  # A temporary hit table (dictionary)
            for line in hit_table:
                h = Hit(sample = sample, hit_line = line, append_sample_name = True)
                if h.query in ht.keys():
                    ht[h.query].append(h)  # >1 hit of the current query sequence is found by BLAST in the current sample
                else:
                    ht[h.query] = [h]
            self.__hit_tables[sample] = ht  # Variable self.__hit_tables[sample_name] is a nested dictionary.
        else:  # No hit is generated at all from the current sample
            self.__hit_counts[sample] = 0
            self.__hit_tables[sample] = None  # A record is created even if the sample does not have any hits, so we won't lose any samples in downstream analysis.
        return

    def count_sample_hits(self, sample):  # Count the number of hits of sample s
        if sample in self.__hit_counts.keys():
            n = self.__hit_counts[sample]
        else:
            print(f"Warning (count_sample_hits): sample {sample} is not found among BLAST outputs", file = sys.stderr)
            n = None
        return n

    def write_hit_sequences(self, query, outdir):
        with open(os.path.join(outdir, '.'.join([query, 'fna'])), 'w') as fasta:  # Override any previous output
            for s, t in self.__hit_tables.items():
                if t != None:
                    if query in t.keys():
                        for h in t[query]:  # Iterate through the list of hits
                            h.write_seq(fasta)
                    else:
                        print(f"Warning (write_hit_sequences): query sequence {query} was not found in sample {s}.", file = sys.stderr)
                else:
                    print(f"Warning (write_hit_sequences): no query sequence was found in sample {s}.", file = sys.stderr)
        return
    
    def compile_tables(self, outdir):
        output_tsv = open(os.path.join(outdir, 'compiled_hits.tsv'), 'w')
        print('\t'.join(['sample'] + HIT_TABLE_COLS), file = output_tsv)  # Print the header line
        for s, t in self.__hit_tables.items():
            if t != None:
                for q, hits in t.items():
                    for h in hits:
                        print('\t'.join([s] + h.attr_values), file = output_tsv)
            else:
                print(f"Warning (compile_tables): no query sequence was found in sample {s}.", file = sys.stderr)
        output_tsv.close()
        return
