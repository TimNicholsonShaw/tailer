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

Local alignment with with FASTA/Q and specific Ensembl IDs of interest or reference fasta
```bash
Tailer -e [comma separated list of EnsIDs no spaces] [FASTA/Q files]
Tailer -f fasta_reference_file.fasta [FASTA/Q files]
```

Optional arguments

* ``-t, --threshold [int, default=100]``
    - Any alignment identified further than this distance in nucleotides from the mature end will be considered spurious and discarded
* ``-x, --trim [int, default=0]``
    - Helper for local mode only, can remove X nucleotides from adapter on the 3' end
* ``-r, --rev_comp``
    - Helper for local mode only. If set, will reverse complement the reads which is necessary for the Lykke-Andersen pipeline
* ``-f, --fasta``
    - Use a fasta file as a reference instead of building one from ensembl IDs (Local Only)
* ``-s, --sequence``
    - Output sequence in the tail file. Useful for debugging.


