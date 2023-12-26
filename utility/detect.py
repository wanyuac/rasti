#!/usr/bin/env python3

"""
MegaBLAST-based search of query sequences against genome assemblies.
Dependencies: BLAST+, Python 3, BioPython, pandas, cd-hit
Example command
    rasti.py --query query/query_genes.fna --genomes *.fna --min_qcov 0 --pause 0.05
Note: this script cannot grep FASTA files for --genomes on Windows OS. Please use Windows's Linux subsystem to run this script.

Copyright (C) 2023-2024 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 14 Jan 2023; the latest update: 27 Dec 2023.
"""

# Standard modules
import os
import sys
import subprocess
from time import sleep

# Project-specific modules
from module.BLAST import BLAST
from module.Queries import Queries
from module.Hit_tables import Hit_tables
from module.SanityCheck import SanityCheck
from module.CDHITEST_cluster_file import CDHITEST_cluster_file

def detect(query, genomes, assembly_suffix, outdir, min_identity, min_qcov, max_evalue, max_match_num, \
           pause, job_reload, cd_hit_est, threads):
    # 1. Environmental settings and sanity checks ###############
    out_dirs = {'root' : outdir, \
                'blast' : os.path.join(outdir, '1_blast'),\
                'parsed' : os.path.join(outdir, '2_parsed'),\
                'extended' : os.path.join(outdir, '3_extended'),\
                'clusters' : os.path.join(outdir, '4_clusters')}
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
    """
    This stage searches query sequences in a multi-FASTA file against each subject genome using megablast/blastn
    and saves each output table (hit_table) as an element in object hit_tables's attribute dictionary
    {target genome : hit table}. Each hit table is a dictionary {query name : one line of hit as stored in a Hit
    object}.
    """
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
    """
    This stage compiles hit tables in object hit_tables into one large table by adding names of subject genomes as keys
    for hit tables. The outcomes of this stage, namely, results of the compile_tables method of hit_tables, consist of
    a TSV file 'compiled_hits.tsv' (in this stage's output directory '2_parsed') and multi-FASTA files '[query name].fna'
    of matched query sequences in subject genomes.
    """
    parsed_out_dir = out_dirs['parsed']
    hit_tables.compile_tables(outdir = parsed_out_dir, extended = False)  # Compile BLAST outputs across all samples into one TSV file
    for q in queries.query_names:  # Create a multi-FASTA file for each query sequence
        hit_tables.write_hit_sequences(query = q, outdir = parsed_out_dir)

    # 4. Extend hits of CDSs to recover alternative start and stop codons ###############
    """
    This stage first walks through each hit in the compiled hit table and recovers alternative start/stop codons in each match
    of CDSs if applicable. No change applies to hits of non-CDS queries. Then the program recompile hit tables and prints the
    compiled hit table into a TSV file 'compiled_hits.tsv' in this stage's output directory '3_extended' and with the same
    column names as those in '2_parsed/compiled_hits.tsv'. In addition, matched query sequences in subject genomes are
    regenerated from updated hit tables and saved in this directory.

    Only an empy file 'no_extended_hit' is created in '3_extended' when no CDS hit is extended.
    """
    ext_out_dir = out_dirs['extended']
    if queries.cds_num > 0:
        hit_tables.extend_cds_hits(subjects = genomes, cds = queries.cds)
        if hit_tables.extension_count > 0:
            hit_tables.compile_tables(outdir = ext_out_dir, extended = True)  # Compile updated hit tables into a large table
            for q in queries.query_names:
                hit_tables.write_hit_sequences(query = q, outdir = ext_out_dir)  # Save updated sequences of hits
            sseq_dir = ext_out_dir  # Input directory (matched subject sequences) of the next stage
        else:
            sseq_dir = parsed_out_dir
    hit_tables.write_extension_records(ext_out_dir)  # Create an empty file 'no_extended_hit' in the output directory if no hit is extended.
    
    # 5. Cluster hits using cd-hit-est ###############
    """ Run cd-hit-est for each FASTA file of matched subject sequences """
    cluster_out_dir = out_dirs['clusters']
    clustr_success = []  # Names of query sequences whose hits were successfully clustered.
    if SanityCheck.cd_hit_est(cd_hit_est):
        for q in queries.query_names:
            f = os.path.join(sseq_dir, '.'.join([q, 'fna']))
            if os.path.exists(f):
                subprocess.run([cd_hit_est, '-i', f, '-o', os.path.join(cluster_out_dir, q + '_representatives.fna'), \
                                '-c', '1.0', '-s', '1.0', '-g', '1', '-d', '0', '-T', threads])
                clustr_file = os.path.join(cluster_out_dir, q + '_representatives.fna.clstr')
                if os.path.exists(clustr_file):
                    clustr_success.append(q)
                    cdhit_output = CDHITEST_cluster_file(clustr_file)
                    cdhit_output.tabulate(os.path.join(cluster_out_dir, q + '_clusters.tsv'))
                else:
                    print(f"Error (sequence clustering): clustering output {clustr_file} was not found.")
            else:
                print(f"Error (sequence clustering): FASTA file of subject sequences matching query {q} was not found. Skipped this file.")
    else:
        print(f"Error (sequence clustering): cd-hit-est could not be found at {cd_hit_est}. No clustering was performed.")

    # Write names of successfully processed query sequences for other sub-commands
    with open(os.path.join(cluster_out_dir, 'clustered_queries.txt'), 'w') as success_record:
        for q in clustr_success:
            success_record.write(q + '\n')
    return
