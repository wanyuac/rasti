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

    def determine_alleles(self, query, compiled_hit_table):
        input_table = os.path.join(self.__input_dir, query + '_clusters.tsv')  # Clustering result of the current query sequence
        if os.path.exists(input_table):
            cluster_table = pandas.read_csv(input_table, sep = '\t')
            cluster_indices = cluster_table['cluster'].unique().tolist()  # A list of integers [0] when there is a single cluster or [0, 1, ...] when there are multiple clusters
            allele_index = 1  # Initial index of non-reference allele
            output_columns = ['query', 'allele'] + cluster_table.columns.tolist()
            allele_assignment = pandas.DataFrame(columns = output_columns)
            for c in cluster_indices:  # Determine if each representative sequence of the clusters is a perfect match of the reference allele of the current gene (sequence cluster)
                table_c = cluster_table[cluster_table['cluster'] == c].copy()  # Take a subset of cluster_table for rows whose cluster == c; output: a data frame of >=1 row; Use the copy method to avoid pandas's SettingWithCopyWarning: https://saturncloud.io/blog/pandas-warning-when-using-map-a-value-is-trying-to-be-set-on-a-copy-of-a-slice-from-a-dataframe/.
                rep_name = str(table_c.loc[table_c['representative'] == 'Y', 'seqid'].iloc[0])  # Extract seqid from the hit corresponding to the cluster's representative sequence
                hit_rep = compiled_hit_table.loc[compiled_hit_table['hit'] == rep_name]  # Extract hit information according to the name of the representative sequence
                hit_rep_similarity = hit_rep['pident'].iloc[0]
                hit_rep_coverage = hit_rep['qcovhsp'].iloc[0]
                if hit_rep_similarity == 100 and hit_rep_coverage == 100:  # A perfect match to the reference allele
                    allele = query  # In this case, use the reference allele's name (gene name) as the allele name, e.g., blaIMP-70
                else:
                    allele = '.'.join([query, str(allele_index)])  # For instance, blaIMP-70.1
                    allele_index += 1
                table_c['allele'] = allele  # Append a column 'allele' to table_c with a constant value
                table_c['query'] = query
                table_c = table_c[output_columns]  # Rearrange columns
                allele_assignment = pandas.concat([allele_assignment, table_c], ignore_index = True)
        else:
            print(f"Error: query-specific cluster table {input_table} does not exist.", file = sys.stderr)
            sys.exit(1)
        return allele_assignment
    
    def create_allele_db(self, query, allele_assignment, outdir):
        """ Rename representative sequences in clusters of query's hits according to the allele assignment from method determine_alleles """
        input_fasta = os.path.join(self.__input_dir, query + '_representatives.fna')
        cluster_representatives = allele_assignment.loc[allele_assignment['representative'] == 'Y']  # The input FASTA file only contains representative sequences.
        if os.path.exists(input_fasta):
            records_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
            output_fasta = open(os.path.join(outdir, query + '_alleles.fna'), 'w')
            for seqid in records_dict.keys():
                allele_info = cluster_representatives.loc[cluster_representatives['seqid'] == seqid]
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

    def update_compiled_hit_table(self, compiled_hit_table, allele_assignments):
        """ Append column 'allele' to the end of the table """
        allele_hit_mapping = allele_assignments[['seqid', 'allele']]
        allele_hit_mapping = allele_hit_mapping.rename(columns = {'seqid' : 'hit'})
        return pandas.merge(compiled_hit_table, allele_hit_mapping, on = 'hit', how = 'left')  # The merge method works in the same way as the merge function in R.
    
    def create_hit_matrix(self, updated_hit_table, genomes, outdir):
        """" Create a genetic matrix with alleles (hits) """
        return
