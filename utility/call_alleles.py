"""
This utility calls alleles of each query sequence (e.g., CDSs) and generates TSV-formatted reports.
It follows the pipeline in script mk_allele_matrix.sh at https://github.com/wanyuac/PAMmaker and
requires Python 3 and the pandas package.

Copyright (C) 2023-2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 27 Dec 2023; the latest update: 15 Apr 2024.
"""

# Standard modules
import os
import sys
import pandas

# Project-specific modules
from module.SanityCheck import SanityCheck
from module.Allele_caller import Allele_caller

def call_alleles(hit_table, representatives_dir, sample_list, outdir):
    SanityCheck.output_dir(outdir)
    if os.path.exists(hit_table) and os.path.exists(sample_list):
        compiled_hit_table = pandas.read_csv(hit_table, sep = '\t')  # Import a compiled hit table, such as compiled_hits_with_extensions.tsv in output subdirectory 3_extended.
        queries = compiled_hit_table['qseqid'].unique().tolist()  # Names of query sequences, namely, reference alleles of query genes / sequence clusters
        if os.path.exists(representatives_dir):
            allele_assignments = dict()
            allele_caller = Allele_caller(indir = representatives_dir, outdir = outdir)
            allele_assignments = pandas.DataFrame(columns = ['query', 'allele', 'cluster', 'index', 'seqid', 'length', 'identity', 'representative'])
            for q in queries:
                q_alleles = allele_caller.determine_alleles(query = q, compiled_hit_table = compiled_hit_table)
                allele_assignments = pandas.concat([allele_assignments, q_alleles], ignore_index = True)  # Concatenate data frames
                allele_caller.create_allele_db(query = q, allele_assignment = q_alleles, outdir = outdir)
                updated_hit_table = allele_caller.update_compiled_hit_table(compiled_hit_table, allele_assignments)
            allele_assignments.to_csv(os.path.join(outdir, 'allele_assignments.tsv'), sep = '\t', index = False)  # Save all allele assignments in a TSV file
            updated_hit_table.to_csv(os.path.join(outdir, 'compiled_hit_table_updated.tsv'), sep = '\t', index = False)
            allele_caller.create_hit_matrix(updated_hit_table, sample_list, sorted(queries), outdir)  # Unlike the list's sort() method, the sorted function returns a new list rather than chaning the original one.
        else:
            print(f"Error: the directory of input FASTA files of representative hit sequences does not exist.", file = sys.stderr)
            sys.exit(1)
    else:
        print(f"Error: the table of hits ({hit_table}) and/or the sample list {sample_list} do not exist.", file = sys.stderr)
        sys.exit(1)
    return
