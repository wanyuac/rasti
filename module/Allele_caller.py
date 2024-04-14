"""
Class Allele_caller for assigning allele identifiers for each gene / sequence cluster.

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 14 Apr 2024; the latest update: 14 Apr 2024.
"""

import os
import sys
import pandas
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Allele_caller:
    """ Assign allele identifiers to a gene or sequence cluster. """
    def __init__(self, indir, outdir):
        self.__input_dir = indir
        self.__outdir = outdir
        return
    
    @property
    def input_dir(self):
        return(self.__input_dir)

    @property
    def output_dir(self):
        return self.__outdir

    def determine_alleles(self, query, compiled_hit_table, outdir):
        input_table = os.path.join(self.__input_dir, query + '_clusters.tsv')
        if os.path.exists(input_table):
            cluster_table = pandas.read_csv(input_table, sep = '\t')
            cluster_indices = cluster_table['cluster'].unique().tolist()  # A list of integers [0] when there is a single cluster or [0, 1, ...] when there are multiple clusters
            allele_index = 1  # Initial index of non-reference allele
            allele_assignment = pandas.DataFrame(columns = ['cluster', 'representative', 'allele'])  # Initiate an empty data frame
            for c in cluster_indices:  # Determine if each representative sequence of the clusters is a perfect match of the reference allele of the current gene (sequence cluster)
                table_c = cluster_table[cluster_table['cluster'] == c]  # Take a subset of cluster_table for rows whose cluster == c; output: a data frame of >=1 row
                rep_name = str(table_c.loc[table_c['representative'] == 'Y', 'seqid'].iloc[0])  # Extract seqid from the hit corresponding to the cluster's representative sequence
                hit_rep = compiled_hit_table.loc[compiled_hit_table['hit'] == rep_name]  # Extract hit information according to the name of the representative sequence
                hit_rep_similarity = hit_rep['pident'].iloc[0]
                hit_rep_coverage = hit_rep['qcovhsp'].iloc[0]
                if hit_rep_similarity == 100 and hit_rep_coverage == 100:  # A perfect match to the reference allele
                    allele = query  # In this case, use the reference allele's name (gene name) as the allele name, e.g., blaIMP-70
                else:
                    allele = '.'.join([query, str(allele_index)])  # For instance, blaIMP-70.1
                    allele_index += 1
                allele_assignment.loc[len(allele_assignment.index)] = [c, rep_name, allele]  # Append a row to the data frame
            allele_assignment.to_csv(os.path.join(outdir, query + '_allele_assignment.tsv'), sep = '\t')
        else:
            print(f"Error: query-specific cluster table {input_table} does not exist.", file = sys.stderr)
            sys.exit(1)
        return allele_assignment
    
    def create_allele_db(self, query, allele_assignment, outdir):
        """ Rename representative sequences in clusters of query's hits according to the allele assignment from method determine_alleles """
        input_fasta = os.path.join(self.__input_dir, query + '_representatives.fna')
        if os.path.exists(input_fasta):
            records_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
            output_fasta = open(os.path.join(outdir, query + '_alleles.fna'), 'w')
            for seqid in records_dict.keys():
                allele_info = allele_assignment.loc[allele_assignment['representative'] == seqid]
                new_seq_name = allele_info['allele'].iloc[0]
                c = allele_info['cluster'].iloc[0]
                seq_record = records_dict[seqid]
                new_seq_descr = seq_record.description[len(seq_record.id) + 1 : ]  # Sequence description is the header line, so this command removes seq_record.id and the following white space from the header.
                new_seq_descr = f'{new_seq_descr}|{seqid}|cd-hit-est_cluster:{c}'
                output_fasta.write(f'>{new_seq_name} {new_seq_descr}\n')
                output_fasta.write(str(seq_record.seq) + '\n')
            output_fasta.close()
        else:
            print(f"Error: FASTA file {input_fasta} does not exist.", file = sys.stderr)
            sys.exit(1)
        return

    def update_compiled_hit_table(self, compiled_hit_table, allele_assignments, outdir):
        pass
        return