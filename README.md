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

Required Arguments

* ``-a, --annotation``
    - GTF formatted database to be used to infer 3' ends. Ensembl GTFs have worked well for this pipeline
* ``[SAM or BAM Files]``
    - Space separated list of sam/bam formatted file locations to perform 3' end tail analysis

Local alignment with with FASTA/Q and specific Ensembl IDs of interest or reference fasta
```bash
Tailer -e [comma separated list of EnsIDs no spaces] [FASTA/Q files]
Tailer -f fasta_reference_file.fasta [FASTA/Q files]
```

Required Arguments

* ``-e, --ensids``
    - Comma separated list of ensembl IDs. Either -e or -f inputs are required.
* ``-f, --fasta``
    - FASTA formatted file to align reads against. Either -f or -e inputs are required
* ``[Trimmed FASTQ files]``
    - Space separated FASTQ file locations for analysis.


Optional arguments

* ``-t, --threshold [int, default=100]``
    - Any alignment identified further than this distance in nucleotides from the mature end will be considered spurious and discarded
* ``-x, --trim [int, default=0]``
    - Helper for local mode only, can remove X nucleotides from adapter on the 3' end
* ``-read, --read [int, default=1]``
    - Paired end only. 1 or 2 to signify which end contains the 3' end information.
* ``-r, --rev_comp``
    - If set, will reverse complement the reads. Necessary for some library prep methods
* ``-f, --fasta``
    - Use a fasta file as a reference instead of building one from ensembl IDs (Local Only)
* ``-s, --sequence``
    - Output sequence in the tail file. Useful for debugging.


