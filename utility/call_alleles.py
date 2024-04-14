"""
This utility calls alleles of each query sequence (e.g., CDSs) and generates TSV-formatted reports.
It follows the pipeline in script mk_allele_matrix.sh at https://github.com/wanyuac/PAMmaker and
requires Python 3 and the pandas package.

Copyright (C) 2023-2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 27 Dec 2023; the latest update: 14 Apr 2024.
"""

# Standard modules
import os
import sys
import pandas

# Project-specific modules
from module.SanityCheck import SanityCheck
from module.Allele_caller import Allele_caller

def call_alleles(hit_table, representatives_dir, outdir):
    SanityCheck.output_dir(outdir)
    if os.path.exists(hit_table):
        compiled_hit_table = pandas.read_csv(hit_table, sep = '\t')  # Import a compiled hit table, such as compiled_hits_with_extensions.tsv in output subdirectory 3_extended.
        queries = compiled_hit_table['qseqid'].unique().tolist()  # Names of query sequences, namely, reference alleles of query genes / sequence clusters
        if os.path.exists(representatives_dir):
            allele_assignments = dict()
            allele_caller = Allele_caller(indir = representatives_dir, outdir = outdir)
            for q in queries:
                alleles = allele_caller.determine_alleles(query = q, compiled_hit_table = compiled_hit_table, outdir = outdir)
                allele_assignments[q] = alleles
                allele_caller.create_allele_db(query = q, allele_assignment = alleles, outdir = outdir)
            updated_hit_table = allele_caller.update_compiled_hit_table(compiled_hit_table, allele_assignments, outdir)
        else:
            print(f"Error: the directory of input FASTA files of representative hit sequences does not exist.", file = sys.stderr)
            sys.exit(1)
    else:
        print(f"Error: the table of hits ({hit_table}) does not exist.", file = sys.stderr)
        sys.exit(1)
    return
