FX Tools
========

This repository contains some tools for processing FASTA and FASTQ files.  They
are implemented using the SeqAn library.

The current status is **experimental**.

fx_convert
----------

Conversion from FASTQ to FASTA and between different quality types of FASTQ
files.

fx_faidx
--------

Indexing of FASTA file and fast indexed access to FASTA files.

This is the equivalent to ``samtools faidx``.

fx_sak
------

Slicing and dicing of FASTA and FASTQ files: Extract certain sequences
by index or name prefix and infixes thereof).  Also allows the
conversion from FASTQ to FASTA, padding of FASTA to FASTQ with dummy
qualities.
