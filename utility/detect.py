#!/usr/bin/env python3

"""
MegaBLAST-based search of query sequences against genome assemblies.
Dependencies: BLAST+, Python 3, BioPython, pandas, cd-hit
Example command
    rasti.py --query query/query_genes.fna --genomes *.fna --min_qcov 0 --pause 0.05
Note: this script cannot grep FASTA files for --genomes on Windows OS. Please use Windows's Linux subsystem to run this script.

Copyright (C) 2023-2024 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 14 Jan 2023; the latest update: 16 Dec 2023.
"""

# Standard modules
import os
import sys
from time import sleep

# Project-specific modules
from module.BLAST import BLAST
from module.Queries import Queries
from module.Hit_tables import Hit_tables
from module.SanityCheck import SanityCheck

def detect(query, genomes, assembly_suffix, outdir, min_identity, min_qcov, max_evalue, max_match_num, pause, job_reload):
    # 1. Environmental settings and sanity checks ###############
    out_dirs = {'root' : outdir, \
                'blast' : os.path.join(outdir, '1_blast'), \
                'parsed' : os.path.join(outdir, '2_parsed'),\
                'extended' : os.path.join(outdir, '3_extended')}
    for d in out_dirs.values():
        SanityCheck.output_dir(d)
    if SanityCheck.output_file(f = query, message = True):
        blast = BLAST(query_fasta = query, \
                      min_identity = SanityCheck.parameter_range(v = min_identity, v_min = 70.0, v_max = 100.0, v_reset = 80.0, n = 'min_identity'),\
                      min_qcov = SanityCheck.parameter_range(v = min_qcov, v_min = 0, v_max = 100.0, v_reset = 80.0, n = 'min_qcov'), max_evalue = max_evalue,\
                      max_hits = SanityCheck.parameter_range(v = max_match_num, v_min = 1, v_max = 200, v_reset = 5, n = 'max_match_num'))
        blast.check_executives()
    else:
        sys.exit(1)
    genomes = SanityCheck.fasta_files(fs = genomes, suf = '.' + assembly_suffix)  # Return: {genome name : Path of the assembly's FASTA file}
    n = len(genomes)  # Number of subject FASTA files
    if n > 0:
        print(f"Number of genomes: {n}", file = sys.stdout)
    else:
        print("Error: none of input genomes exists.", file = sys.stderr)
        sys.exit(1)
    delay_sec = pause  # Number of seconds to hold between consecutive BLAST searches
    delay_iterations = (delay_sec > 0 and delay_sec < 60)

    # 2. Iteratively run megaBLAST through genomes ###############
    blast_out_dir = out_dirs['blast']
    queries = Queries(fasta = query)
    queries.write_query_lengths(tsv = os.path.join(out_dirs['root'], 'query_lengths.tsv'))  # Save a two-column table of lengths of query sequences
    hit_tables = Hit_tables()  # Initialise a Hit_tables object
    if job_reload:
        for g in genomes.keys():
            hit_tables.add_table(sample = g, hit_table = blast.read(subject_name = g, input_dir = blast_out_dir))  # Import existing BLAST outputs
    elif delay_iterations:  # Do fresh BLAST searches
        for g, fasta in genomes.items():
            hit_tables.add_table(sample = g, hit_table = blast.search(subject_name = g, subject_fasta = fasta, outdir = blast_out_dir))  # Method blast.search may return None when no hit is found in a subject genome.
            sleep(delay_sec)
    else:
        for g, fasta in genomes.items():
            hit_tables.add_table(sample = g, hit_table = blast.search(subject_name = g, subject_fasta = fasta, outdir = blast_out_dir))
    
    # 3. Parse and compile BLAST results ###############
    parsed_out_dir = out_dirs['parsed']
    hit_tables.compile_tables(outdir = parsed_out_dir, extended = False)  # Compile BLAST outputs across all samples into one TSV file
    for q in queries.query_names:  # Create a multi-FASTA file for each query sequence
        hit_tables.write_hit_sequences(query = q, outdir = parsed_out_dir)

    # 4. Extend hits of CDSs to recover alternative start and stop codons ###############
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
