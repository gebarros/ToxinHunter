# ToxinHunter

This repository contains a bioinformatics pipeline to identify potential snake venomon toxins from de novo transcriptome assembly.

### Requirements:
    - Python3
    - Blast+
    - TransDecoder
    - BioPython

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

### * [Specifications](specifications.md)
