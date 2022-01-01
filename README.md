# BioSeq

Contains the BioSeq class which will perform some standard bioinformatics processing tools. 

## BioSeq
Functionality includes:

- Provides the information of the sequence including the sequence, sequence type, label, and length
- Generates a random sequence with a provided length. Can create DNA or RNA sequences.
- Counts nucleotide frequency and returns a disctionary
- Performs DNA -> RNA transcription
- Creates the reverse complement of the provided sequence
- Calculates the GC content in a sequence 
- Calculates the GC content in subsests of sequences
- Translates a sequence
- Calculates the frequence of each codons encoding for a given amino acid in a sequence
- Generates six reading frames and computes possible proteins from reading frames
- Computes all possible proteins for all open reading frames
- Calculates hamming distance of two provide sequences

## Utilities
Functionality includes:
- Will read and write to a text file. 
- Will open and read FASTA files

## BioStructs

- Contains nucleotide bases
- DNA codons
- RNA codons

## main
- An example main file
