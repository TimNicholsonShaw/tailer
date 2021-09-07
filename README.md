# Tailer: A tool for 3' end analysis of non-polyadenylated RNAs from 3' sequencing data

## Requirements
- ``Python >= 3.6``
- ``pysam``
- ``Bio``
- ``gffuitls``
- ``requests``
- ``tqdm``

## Installation

```bash
pip install jla-tailer
```

## Usage

Global alignment with a GTF annotation and SAM/BAM formatted files
```bash
Tailer -a [GTF Annotation] [SAM or BAM Files]
```

Local alignment with with FASTA/Q and specific Ensembl IDs of interest
```bash
Tailer -e [comma separated list of EnsIDs] [FASTA/Q files]
```

Optional arguments

* ``-t, --threshold [int, default=100]``
    - Any identified further than this distance in nucleotides from the mature end will be considered spurious and discarded


