#!/usr/bin/env python

"""
This module defines class Hits.

Dependencies: Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 25 Dec 2023.
"""

import os
import sys
import subprocess
import pandas as pd
from copy import deepcopy
from module.Hit import Hit, HIT_ATTRS

class Hit_tables:
    """
    Manage raw BLAST outputs as a dictionary {genome : BLAST's output table, either with CDS hits
    extended to recover alternative start/stop codons}
    """
    def __init__(self):
        self.__hit_tables = dict()  # A nested dictionary {genome : table} of raw BLAST outputs, where 'table' is a dictionary of lists of Hit objects
        self.__extensions = pd.DataFrame(columns = ['Hit', 'Genome', 'Action'])
        return
    
    @property
    def table_num(self):  # Number of hit tables
        return len(self.__hit_tables)

    @property
    def sample_names(self):
        return list(self.__hit_tables.keys())
    
    @property
    def extensions(self):
        return self.__extensions
    
    @property
    def extension_num(self):
        return self.__extensions.shape[0]  # Number of rows

    def add_table(self, sample, hit_table):  # Parameter hit_table is a list of rows in BLAST's raw tabulated output.
        if hit_table != None:  # Each element of the dictionary __hit_tables stores information of rows in the raw output TSV file when the output is not empty. Otherwise, the table is None.
            ht = dict()  # A temporary hit table (dictionary)
            hits_num = dict()
            for line in hit_table:
                h = Hit(sample = sample, hit_line = line, append_sample_name = True)
                q = h.query
                if q in ht.keys():  # More than one hit of the current query sequence is found by BLAST in the current sample (hit table)
                    n = hits_num[q]
                    if n == 1:  # There is only a single hit in list ht[h.query] by far.
                        h_prev = ht[q][0]
                        h_prev.id = ':'.join([h_prev.id, '1'])  # For example, 'gene1@sample1:1'
                        ht[q] = [h_prev]
                    n += 1
                    h.id = ':'.join([h.id, str(n)])  # For example, 'gene1@sample1:2' for the second hit of gene 1 in sample 1
                    ht[q].append(h)
                    hits_num[q] = n
                else:
                    ht[q] = [h]
                    hits_num[q] = 1
            self.__hit_tables[sample] = ht  # Variable self.__hit_tables[sample_name] is a nested dictionary.
        else:  # No hit is generated at all from the current sample
            self.__hit_tables[sample] = None  # A record is created even if the sample does not have any hits, so we won't lose any samples in downstream analysis.
        return

    def write_hit_sequences(self, query, outdir):
        with open(os.path.join(outdir, '.'.join([query, 'fna'])), 'w') as fasta:  # Override any previous output
            for s, t in self.__hit_tables.items():  # Iterate through hit tables by sample names
                if t != None:
                    if query in t.keys():
                        for h in t[query]:  # Iterate through the list of hits
                            h.write_seq(fasta)
                    else:
                        print(f"Warning (write_hit_sequences): query sequence {query} was not found in sample {s}.", file = sys.stderr)
                else:
                    print(f"Warning (write_hit_sequences): no query sequence was found in sample {s}.", file = sys.stderr)
        return
    
    def compile_tables(self, outdir, extended):
        """ Concatenate raw output tables of megaBLAST """
        attrs = deepcopy(HIT_ATTRS[0 : (len(HIT_ATTRS) - 2)]) if extended else deepcopy(HIT_ATTRS)  # Remove 'evalue' and 'bitscore' when extended = True; Deepcopy must be used. Otherwise, HIT_ATTRS will be changed by the append method (https://stackoverflow.com/questions/24345712/python-list-of-objects-changes-when-the-object-that-was-input-in-the-append-f)
        attrs.append('hslen')  # Use attrs = attrs + ['hslen'] if deepcopy isn't used.
        output_tsv = open(os.path.join(outdir, 'compiled_hits.tsv'), 'w')
        print('\t'.join(['sample', 'hit', 'qseqid', 'sseqid'] + attrs), file = output_tsv)  # Print the header line
        for s, t in self.__hit_tables.items():  # Iterate through hit tables by sample names
            if t != None:
                for hits in t.values():
                    for h in hits:
                        print('\t'.join([s, h.id, h.query, h.contig] + h.attr_values(extended)), file = output_tsv)
            else:
                print(f"Warning (compile_tables): no query sequence was found in sample {s}.", file = sys.stderr)
        output_tsv.close()
        return

    def extend_hits(self, subjects, cds):
        """
        Extend CDSs to recover alternative start and stop codons
        Parameters:
          subjects, a dictionary {genome name : path to the FASTA file}, which is generated by function check_assemblies
          cds, a list of CDS names from a Queries object
        """
        extension_table = []
        for genome, hit_table in self.__hit_tables.items():  # Go through every hit table
            for query, hit_list in hit_table.items():  # Go through every hit in the current hit table
                if query in cds:
                    new_hit_list = list()
                    for hit in hit_list:
                        extension_action = hit.extend_cds(f = subjects[genome])
                        if extension_action != "":
                            extension_table.append([hit.id, genome, extension_action])
                        new_hit_list.append(hit)
                    self.__hit_tables[genome][query] = new_hit_list
        self.__extensions = pd.DataFrame(data = extension_table, columns = ['Hit', 'Genome', 'Action'])
        return

    def write_extension_records(self, outdir):
        if self.__extensions.shape[0] > 0:
            tsv = os.path.join(outdir, 'extension_records.tsv')
            self.__extensions.to_csv(tsv, sep = '\t', index = False)
        else:
            tsv = os.path.join(outdir, 'no_extended_hit')
            subprocess.run(f'touch {tsv}', shell = True)
        return
