#!/usr/bin/env python

"""
This module defines class Hit_tables, which parses BLAST's outputs across samples and store them as individual tables of hits.

Dependencies: Python 3, pandas

Copyright (C) 2023 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 31 Jan 2025.
"""

import os
import sys
import subprocess
import pandas
from copy import deepcopy
from module.Hit import Hit, HIT_ATTRS

class Hit_tables:
    """
    Manage raw BLAST outputs as a dictionary of tables:
    {sample : BLAST's output table, either with CDS hits extended to recover alternative start/stop codons}
    """
    def __init__(self):
        self.__hit_tables = dict()  # A nested dictionary {sample : table} of raw BLAST outputs, where 'table' is a dictionary of lists of Hit objects
        self.__extensions = pandas.DataFrame(columns = ['Hit', 'Sample', 'Action'])
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
    def extension_count(self):  # Number of CDSs extended
        return self.__extensions.shape[0]  # Number of rows

    def add_table(self, sample, hit_table):
        """
        This method parses megablast/blastn's tab-delimited output table of one sample (hit_table) and saves
        it as a dictionary object {query name : hit information as stored as a Hit object}.
        
        This method deals with multiple hits of the same query sequence. Specifically, when >1 hits are present
        for a query sequence q, the query name follows the format '[gene name]@[sample name]:[index]'; otherwise,
        the query name follows '[gene name]@[sample name]'.
        
        Parameter hit_table of this method is a list of rows in BLAST's raw tabulated output. This list is the output of the serach method
        of class BLAST.
        """
        if hit_table != None:  # Each element of the dictionary __hit_tables stores information of rows in the raw output TSV file when the output is not empty. Otherwise, the table is None.
            ht = dict()  # A temporary hit table (dictionary)
            hits_num = dict()
            for line in hit_table:
                h = Hit(sample = sample, hit_line = line, append_sample_name = True)
                q = h.query
                if q in ht.keys():  # More than one hit of the current query sequence is found by BLAST in the current sample (hit table)
                    n = hits_num[q]
                    if n == 1:  # There is only a single hit in list ht[h.query] by far.
                        h_prev = ht[q][0]  # An Hit object
                        h_prev.id = ':'.join([h_prev.id, '1'])  # For example, 'gene1@sample1:1'
                        ht[q] = [h_prev]  # Override ht[q][0] with an Hit object of an updated query name 'gene1@sample1:1'
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
    
    def compile_tables(self, outdir, extended):
        """ Concatenates raw output tables of megaBLAST and prints results into a TSV file """
        if extended:
            """
            Deepcopy must be used in the following commands. Otherwise, HIT_ATTRS will be changed by the append method:
            https://stackoverflow.com/questions/24345712/python-list-of-objects-changes-when-the-object-that-was-input-in-the-append-f.
            """
            attrs = deepcopy(HIT_ATTRS[0 : (len(HIT_ATTRS) - 2)])  # Remove 'evalue' and 'bitscore' when extended = True
            output_tsv = open(os.path.join(outdir, 'compiled_hits_with_extensions.tsv'), 'w')
        else:
            attrs = deepcopy(HIT_ATTRS)
            output_tsv = open(os.path.join(outdir, 'compiled_hits.tsv'), 'w')
        attrs.append('hslen')  # Length of the original hit. hslen = send - sstart + 1 according to nucleotide BLAST's output. 
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

    def write_hit_sequences(self, query, outdir):
        """ This method creates a multi-FASTA file of matched sequences of the query in subject samples. """
        with open(os.path.join(outdir, '.'.join([query, 'fna'])), 'w') as fasta:  # Override any previous output
            for s, t in self.__hit_tables.items():  # Iterate through hit tables by sample names
                if t != None:
                    if query in t.keys():
                        for h in t[query]:  # Iterate through the list of hits
                            h.write_seq(fasta)
                else:
                    print(f"Warning (write_hit_sequences): no query sequence was found in sample {s}.", file = sys.stderr)
        return
    
    def extend_cds_hits(self, subjects, cds):
        """
        This method extends CDSs (specified by parameter cds) in each hit table, if any, to recover alternative start
        and stop codons. No change applies to other queries in this table.

        Parameters:
          - subjects, a dictionary {sample name : path to the FASTA file}, which is generated by function check_assemblies
          - cds, a list of CDS names in a Queries object's attribute 'cds'
        """
        extension_table = []
        for sample, hit_table in self.__hit_tables.items():  # Go through every hit table, namely, BLAST outputs of each sample
            if hit_table != None:
                for query, hit_list in hit_table.items():  # Go through every hit in the current hit table
                    if query in cds:
                        new_hit_list = list()
                        for hit in hit_list:
                            extension_action = hit.extend_cds(f = subjects[sample])
                            if extension_action != "":
                                extension_table.append([hit.id, sample, extension_action])
                            new_hit_list.append(hit)
                        self.__hit_tables[sample][query] = new_hit_list  # Modify the corresponding hit table; no change to other tables that do not match these if conditions
        """ Create a simple table for an overview of all extended CDSs. This information is also added to compiled_hits_extended.tsv. """
        self.__extensions = pandas.DataFrame(data = extension_table, columns = ['Hit', 'Sample', 'Action'])
        return

    def write_extension_records(self, outdir):
        if self.__extensions.shape[0] > 0:
            tsv = os.path.join(outdir, 'extension_records.tsv')
            self.__extensions.to_csv(tsv, sep = '\t', index = False)
        else:
            tsv = os.path.join(outdir, 'no_extended_hit')
            subprocess.run(f'touch {tsv}', shell = True)
        return
