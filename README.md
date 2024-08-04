# Iterative assembly-based search for target nucleotides

<img src="logo/rasti.png" alt="rasti logo" style="float: left; margin-right: 10px; width: 20%;" />

Current release: `rasti v0.0.1`
Latest documentation update: 4 August 2024



## 1. Overview

Rasti (**i**terative **a**ssembly-based **s**earch for **t**arget nucleot**i**des) is a Python wrapper of nucleotide BLAST that implements three sequential utilities to screen assemblies for any bacterial/archaeal/plant plastid genes:  

* `detect`: searching query sequences against contigs in a genome/metagenome assembly — here, each query is the reference allele of a gene or locus;  
* `call_alleles`: taking as input outcomes of the `detect` method, this utility assigns allele identifiers to hits of each query sequence and generates a matrix of allelic presence-absence across samples;  
* `aln2mut` (alignment to mutations): identifying mutations from the global alignment of a gene's alleles.  

Before using rasti, users are recommended to reduce query redundancy by assigning similar queries (for example, using [CD-HIT-EST](https://github.com/weizhongli/cdhit) to cluster sequences based on a minimum nucleotide identity and coverage of 80%) into clusters (often refer to as genes) and selecting representative sequences of these clusters for sequence search, which is the same as [SRST2](https://github.com/katholt/srst2) and [ARIBA](https://github.com/sanger-pathogens/ariba) do.



### How rasti works

Users are expected to run rasti's `detect`, `call_alleles`, and `aln2mut` sequentially.

* `detect`: this method searches query sequences (from `queries.fna` for example) against each sample's assembly (namely, the subject) using megaBLAST and reads the result table. Then it compiles result tables across all samples into a single table. Since megaBLAST usually cannot identify the complete coding sequence (CDS) when alternative start/stop codons are present in the subject sequence compared with the query allele, rasti attempts to overcome this technical limitation by expanding each hit of a CDS in the subject sequence and searching for any alternative start/stop codons following the [Translation Codon Table 11](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG11) — for this reason, rasti only works for bacterial/archaeal/plant plastid genes but in theory, it would enable users to select the codon table for other organisms. Finally, it clusters hits of each query sequence (each "gene" or "locus") using `cd-hit-est` to identify unique alleles under 100% nucleotide identity and query coverage.
* `call_alleles`: this method assigns a unique identifier to each allele of a gene and creates a allelic matrix that is similar to the outcome of SRST2 and ARIBA. The allele identifier is the same as the query's name if the allele is identical to the reference, and an arbitrary index 1, 2, 3, ... is added to the query's name otherwise. For example, a precise match to the reference allele of gene blaIMP-70 has an allele identifier of blaIMP-70, and identifier blaIMP-70.1 is given to an allele that differs from the reference. The method also pools all alleles of each gene (including the reference allele regardless whether it is present in any sample for the convenience of downstream mutation identification and protein-level comparisons) into one FASTA file with the name `[gene name]_alleles.fna` (*e.g.,* `blaIMP-70_alleles.fna`).
* `aln2mut`: to use this method, users need to run their preferred alignment tool (*e.g.,* MUSCLE or Clustal Omega) to generate a global alignment of alleles for each gene of interest. Then the alignment is saved as a multi-FASTA file and used as input for `aln2mut`, which will report mutations in a VCF file, a matrix, and optionally, a list. The method may also produce an alignment of only variable sites (similar to Sanger Institute's tool [snp-sites](https://github.com/sanger-pathogens/snp-sites)).
  
  

### Dependencies

Rasti does not strongly rely on particular versions of dependencies. Here, I list versions that have been tested for `rasti` during its development.

* [Python 3](https://www.python.org/downloads/) (v3.9.19)  
* [BioPython](https://github.com/biopython/biopython) (v1.83)  
* [Pandas](https://pandas.pydata.org/) (v2.2.2)  
* [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (v2.15.0)  
* [CD-HIT](https://github.com/weizhongli/cdhit) (v4.8.1)  
  
  

## 2. Installation

[Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and [mamba ](https://github.com/mamba-org/mamba) are recommended ways for installing rasti.

```bash
conda create -n rasti python=3.9
conda activate rasti
conda install -c conda-forge biopython
conda install -c conda-forge pandas
conda install -c conda-forge -c bioconda cd-hit
conda install -c conda-forge -c bioconda blast
git clone https://github.com/wanyuac/rasti.git  # cd to your preferred directory first
```



## 3. Quick start

This example assumes `queries.fna` has been prepared (see the next section about preparation of query FASTA files) and a rasti's conda environment has been activated.

```bash
cd ~/analysis  # Or any other working directory

~/bin/rasti/rasti.py detect --query input/queries.fna --assemblies input/ERR*.fna --assembly_suffix fna --outdir output/test --min_identity 90 --min_qcov 90 --cd_hit_est ~/anaconda3/envs/rasti/bin/cd-hit-est --threads 8

~/bin/rasti/rasti.py call_alleles --compiled_hit_table output/test/3_extended/compiled_hits_with_extensions.tsv --sample_list output/test/sample_list.txt --queries input/queries.fna --representatives_dir output/test/4_clusters --outdir output/test/5_alleles

mkdir output/test/6_mutations
~/bin/rasti/rasti.py aln2mut --input output/test/5_alleles/gene1_alleles.fna --outdir output/test/6_mutations --output_prefix gene1_mutations --ref_name gene1 --list --var
```



## 4. Preparing inputs

Rasti takes as input two types of FASTA files: one type for assemblies and the other type for query sequences.



### FASTA files of assemblies

Rasti assumes filenames of assemblies' FASTA files follow the format `[sample name].[fna/fasta]`. For example, `sample1.fna` is a valid filename.



### A FASTA file of query sequences

This input is a multi-FASTA file with a header format in which each coding sequence (CDS) is indicated by keyword "CDS" at the beginning of the annotation field in the sequence header. For example, such sequence headers can be `>seq1 CDS` and `>seq2 CDS|annotations`, and so forth. Hit will not be extended to correct partial matches of alternative start/stop codons if the query is not a CDS, hence an empty output subdirectory `3_extended`.



## 5. `rasti.py detect` method

### Parameters

* `--queries` / `-q`: mandatory input: a multi-FASTA file of query DNA sequences. For coding sequences, add CDS to the beginning of sequence annotations and separated from other annotations with a '|' character in this FASTA file. For example, '>seq CDS|other annotations';

* `--assemblies` / `-a`: mandatory input: FASTA files of assemblies against which queries will be searched;

* `--assembly_suffix` / `-s`: filename extension (fasta/fna/fa, etc) to be removed from assembly filenames in order to get a sample name (default: "fna");

* `--outdir` / `-o`: output directory (default: "output");

* `--min_identity` / `-mi`: minimum percent nucleotide identity for BLAST to identify a match (default: 90.0; range: 70-100);

* `--min_qcov` / `-mq`: minimum percent query coverage for BLAST to identify a match (default: 90.0; range: 0-100);

* `--max_evalue` / `-me`: maximum E-value for BLAST to identify a match (default: 1e-5);

* `--max_match_num` / `-mh`: maximum number of matches reported by BLAST for each query sequence (default: 5; Range: 1-500);

* `--pause` / `-p`: number of seconds to hold between consecutive BLAST searches (default: 0.2; range: 0-60);

* `--reload` / `-r`: flag this option to enable importing existing BLAST outputs without reruning the BLAST search — the pause is disabled in this case;

* `--cd_hit_est` / `-c`: full path of program cd-hit-est;

* `--threads` / `-t`: number of threads for BLAST and CD-HIT-EST.
  
  

## 6. `rasti.py call_alleles` method

### Parameters

* `--compiled_hit_table` / `-t`: compiled hits in the outputs of method `detect`;

* `--sample_list` / `-s`: a test file listing names of sample. It can be "sample_list.txt" in the output directory of `rasti detect;

* `--queries` / `-q`: the same multi-FASTA file of query DNA sequences used for the `detect` method;

* `--representatives_dir` / `-r`: directory of input FASTA files of representatives allele sequences (`*_representatives.fna`);

* `--outdir` / `-o`: output directory (default: "output/5_alleles").
  
  

## 7. `rasti.py aln2mut` method

This method works on both nucleotide and amino acid sequence alignments. Because of synonymous mutations, it is desirable to cluster identical translated sequences (namely, protein sequences of 100% amino acid identity and coverage) using [`CD-HIT`](https://github.com/weizhongli/cdhit) to remove alignment redundancy (namely, run `rasti.py aln2mut` only through representative amino acid sequences from CD-HIT) and ensure the resulting amino acid alignment always shows variation.

### Parameters

* `--input` / `-i`: input alignment in the FASTA format;

* `--outdir` / `-o`: output directory (default: current working directory);

* `--output_prefix` / `-p`: prefix for output files (default: "mutations");

* `--ref_name` / `-r`: name of the reference sequence in the alignment;

* `--list` / `-l`: create a list of alterations in a conventional format (*e.g.*, W25N);

* `--var` / `-v`: create a FASTA-format alignment file of variable sites only.
  
  

## 8. Etymology

"Rasti" is a Lithuanian verb and noun meaning "(to) find" and "(to) discover".



## 9. Development history

Rasti is a combination and enhancement of [NITREc](https://github.com/wanyuac/NITREc/tree/master/Script), [geneDetector](https://github.com/wanyuac/geneDetector), [PAMmaker](https://github.com/wanyuac/PAMmaker), and [aln2mut](https://github.com/wanyuac/aln2mut).
