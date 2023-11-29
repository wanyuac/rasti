# Rasti
Release: TBC
Latest update: 29/11/2023

**R**ecursive **a**ssembly-based **s**earch for **t**arget nucleot**i**des  

This software is a combination and enhancement of [NITREc](https://github.com/wanyuac/NITREc/tree/master/Script) and [geneDetector](https://github.com/wanyuac/geneDetector).  

## Installation
```bash
conda create -n rasti python=3.9
conda activate rasti
conda install -c bioconda blast
conda install -c conda-forge biopython
conda install -c anaconda pandas
conda install -c bioconda cd-hit
git clone https://github.com/wanyuac/rasti.git
```

### Dependencies
Recommend installing latest versions of the following dependencies. I noted versions that have been tested for `rasti` although other versions may also work.
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (v2.15.0)
- [Python 3](https://www.python.org/downloads/) (v3.12.0)
- [BioPython](https://github.com/biopython/biopython) (v1.81)
- [Python module pandas](https://pandas.pydata.org/) (v2.1.3)
- [CD-HIT](https://github.com/weizhongli/cdhit) (v4.8.1)

## Parameters
```bash
python rasti/rasti.py --help

usage: rasti.py [-h] --query QUERY --genomes GENOMES [GENOMES ...] [--assembly_suffix ASSEMBLY_SUFFIX] [--outdir OUTDIR] [--min_identity MIN_IDENTITY] [--min_qcov MIN_QCOV]
                [--max_evalue MAX_EVALUE] [--max_match_num MAX_MATCH_NUM] [--pause PAUSE] [--reload]

Targeted gene detection for assemblies

options:
  -h, --help            show this help message and exit
  --query QUERY, -q QUERY
                        (Mandatory) a multi-Fasta file of query DNA sequences
  --genomes GENOMES [GENOMES ...], -g GENOMES [GENOMES ...]
                        (Mandatory) Fasta files of genome assemblies against which queries will be searched
  --assembly_suffix ASSEMBLY_SUFFIX, -s ASSEMBLY_SUFFIX
                        Filename extension (fasta/fna/fa, etc) to be removed from assembly filenames in order to get a sample name (Default: fna)
  --outdir OUTDIR, -o OUTDIR
                        Output directory (Default: results)
  --min_identity MIN_IDENTITY, -mi MIN_IDENTITY
                        Minimum percent nucleotide identity for BLAST to identify a match (Default: 80.0; range: 70-100)
  --min_qcov MIN_QCOV, -mq MIN_QCOV
                        Minimum percent query coverage for BLAST to identify a match (Default: 80.0; range: 0-100)
  --max_evalue MAX_EVALUE, -me MAX_EVALUE
                        Maximum E-value for BLAST to identify a match (Default: 1e-5)
  --max_match_num MAX_MATCH_NUM, -mh MAX_MATCH_NUM
                        Maximum number of matches reported by BLAST for each query sequence (Default: 5; Range: 1-500)
  --pause PAUSE, -p PAUSE
                        Seconds to be paused between BLAST searches (Default: 0.2; range: 0-60)
  --reload, -r          Flag this option to enable importing existing BLAST outputs without reruning the BLAST search (Option --pause is disabled in this case)
```

## Header format of query sequences
For each coding sequence (CDS):

```fasta
>seq1 CDS

>seq2 CDS|annotations
```
"CDS" is reserved for specifying a CDS. Hit will not be extended if the query is not a CDS, hence an empty output subdirectory `3_extended`.

## Etymology
"Rasti" is a Lithuanian verb and noun meaning "(to) find" and "(to) discover".  
