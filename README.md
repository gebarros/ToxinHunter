# ToxinHunter

This repository contains a bioinformatics pipeline to identify potential snake venomon toxins from de novo transcriptome assembly.

### Requirements:
```
    - Python3
    - Blast+
    - TransDecoder
    - BioPython
```

### Installation using miniconda:

  1- [Install miniconda](https://conda.io/miniconda.html)
  2 - Create conda env using Python3:
  ```
   conda create -n <env_name> python=3
  ```

  3- Install packages using conda:
  ```
   conda install biopython
   conda install -c bioconda transdecoder
   conda install -c bioconda blast
  ```

### Usage:
```
usage: toxinHunter.py [-h] [f] [t] [cov] [id]

Toxin Hunter: find potential toxin from de novo transcriptome assembly

positional arguments:
  f           Fasta file of transcriptome assembly.
  t           Toxins dataset (Fasta format with a *specific header).
  cov         Percentage of coverage used to select contigs in tblastn result
              (Value must be between 0 and 100).
  id          Sample name.

optional arguments:
  -h, --help  show this help message and exit
```

### Output:

 - Potential snake toxins identified from RNAseq assembly (Nucleotide ORFs, Amino acids ORFs and Whole Contig)


[Click here to see specifications.](specifications.md)

 **Contact:** gebarrosbio@gmail.com
