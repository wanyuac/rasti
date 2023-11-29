#!/usr/bin/env python

"""
MegaBLAST-based search of query sequences against genome assemblies.

Dependencies: BLAST+, Python 3, BioPython, pandas, cd-hit

Example command: rasti.py --query query/query_genes.fna --genomes *.fna --min_qcov 0 --pause 0.05

Note: this script cannot grep FASTA files for --genomes on Windows OS. Please use Windows's Linux subsystem to run this script.

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 14 Jan 2023; the latest update: 29 Nov 2023.
"""

import os
import sys
from time import sleep
from argparse import ArgumentParser
from lib.Queries import Queries
from lib.BLAST import BLAST
from lib.Hit_tables import Hit_tables
from lib.utilities import check_dir, check_file, check_assemblies, check_values


def parse_arguments():
    parser = ArgumentParser(description = "Targeted gene detection for assemblies")
    parser.add_argument('--query', '-q', dest = 'query', type = str, required = True, help = "(Mandatory) a multi-Fasta file of query DNA sequences. For coding sequences, add CDS to the beginning of sequence annotations")
    parser.add_argument('--genomes', '-g', nargs = '+', dest = "genomes", type = str, required = True, help = "(Mandatory) Fasta files of genome assemblies against which queries will be searched")
    parser.add_argument('--assembly_suffix', '-s', dest = 'assembly_suffix', type = str, required = False, default = 'fna', help = "Filename extension (fasta/fna/fa, etc) to be removed from assembly filenames in order to get a sample name (Default: fna)")
    parser.add_argument('--outdir', '-o', dest = 'outdir', type =str, required = False, default = 'results', help = "Output directory (Default: results)")
    parser.add_argument('--min_identity', '-mi', dest = 'min_identity', type = float, default= 80.0, required = False, help = "Minimum percent nucleotide identity for BLAST to identify a match (Default: 80.0; range: 70-100)")
    parser.add_argument('--min_qcov', '-mq', dest = 'min_qcov', type = float, default = 80.0, required = False, help = "Minimum percent query coverage for BLAST to identify a match (Default: 80.0; range: 0-100)")
    parser.add_argument('--max_evalue', '-me', dest = 'max_evalue', type = str, default = '1e-5', required = False, help = "Maximum E-value for BLAST to identify a match (Default: 1e-5)")
    parser.add_argument('--max_match_num', '-mh', dest = 'max_match_num', type = int, default = 5, required = False, help = "Maximum number of matches reported by BLAST for each query sequence (Default: 5; Range: 1-500)")
    parser.add_argument('--pause', '-p', dest = 'pause', type = float, default = 0.2, required = False, help = "Seconds to be paused between BLAST searches (Default: 0.2; range: 0-60)")
    parser.add_argument('--reload', '-r', dest = 'reload', action = 'store_true', help = "Flag this option to enable importing existing BLAST outputs without reruning the BLAST search (Option --pause is disabled in this case)")
    return parser.parse_args()


def main ():
    args = parse_arguments()

    # Environmental settings and sanity check
    out_dirs = {'root' : args.outdir, 'blast' : os.path.join(args.outdir, '1_blast'), 'parsed' : os.path.join(args.outdir, '2_parsed'),\
                'extended' : os.path.join(args.outdir, '3_extended')}
    for d in out_dirs.values():
        check_dir(d)
    if check_file(f = args.query, message = True):
        blast = BLAST(query_fasta = args.query, min_identity = check_values(v = args.min_identity, v_min = 70.0, v_max = 100.0, v_reset = 80.0, n = 'min_identity'),\
                      min_qcov = check_values(v = args.min_qcov, v_min = 0, v_max = 100.0, v_reset = 80.0, n = 'min_qcov'), max_evalue = args.max_evalue,\
                      max_hits = check_values(v = args.max_match_num, v_min = 1, v_max = 200, v_reset = 5, n = 'max_match_num'))
        blast.check_executives()
    else:
        sys.exit(1)
    genomes = check_assemblies(fs = args.genomes, suf = '.' + args.assembly_suffix)  # Return: {genome name : Path of the assembly's FASTA file}
    n = len(genomes)
    if n > 0:
        print(f"Number of genomes: {n}", file = sys.stdout)
    else:
        print("Error: none of input genomes exists.", file = sys.stderr)
        sys.exit(1)
    delay_sec = args.pause
    delay_iterations = (delay_sec > 0 and delay_sec < 60)

    # Iteratively run megaBLAST through genomes
    blast_out_dir = out_dirs['blast']
    queries = Queries(fasta = args.query)
    queries.write_query_lengths(tsv = os.path.join(out_dirs['root'], 'query_lengths.tsv'))  # Save a two-column table of lengths of query sequences
    hit_tables = Hit_tables()  # Initialise a Hit_tables object
    if args.reload:
        for g in genomes.keys():
            hit_tables.add_table(sample = g, hit_table = blast.read(subject_name = g, input_dir = blast_out_dir))  # Import existing BLAST outputs
    elif delay_iterations:  # Do fresh BLAST searches
        for g, fasta in genomes.items():
            hit_tables.add_table(sample = g, hit_table = blast.search(subject_name = g, subject_fasta = fasta, outdir = blast_out_dir))  # Method blast.search may return None when no hit is found in a subject genome.
            sleep(delay_sec)
    else:
        for g, fasta in genomes.items():
            hit_tables.add_table(sample = g, hit_table = blast.search(subject_name = g, subject_fasta = fasta, outdir = blast_out_dir))
    
    # Parse and compile BLAST results
    parsed_out_dir = out_dirs['parsed']
    hit_tables.compile_tables(outdir = parsed_out_dir, extended = False)  # Compile BLAST outputs across all samples into one TSV file
    for q in queries.query_names:  # Create a multi-FASTA file for each query sequence
        hit_tables.write_hit_sequences(query = q, outdir = parsed_out_dir)

    # Extend hits of CDSs to recover alternative start and stop codons
    cds = queries.cds
    if len(cds) > 0:
        ext_out_dir = out_dirs['extended']
        hit_tables.extend_hits(subjects = genomes, cds = cds)
        if hit_tables.extension_num > 0:
            hit_tables.compile_tables(outdir = ext_out_dir, extended = True)
        hit_tables.write_extension_records(ext_out_dir)  # Creates an empty file 'no_extended_hit' in the output directory if no hit is extended.
        for q in queries.query_names:
            hit_tables.write_hit_sequences(query = q, outdir = ext_out_dir)
    return


if __name__ == '__main__':
    main()
