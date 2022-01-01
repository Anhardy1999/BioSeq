# Toolset/Code testing file

from os import write
from BioSeq  import BioSeq
from utilities import read_FASTA, readTextFile, writeTextFile

test_dna = BioSeq()
test_dna.generate_rand_seq(40, 'DNA')
print(test_dna.get_seq_info())
print(test_dna.nucleotide_freq())
print(test_dna.transcription())
print(test_dna.reverse_complement())
print(test_dna.gc_content())
print(test_dna.gc_content_subset())
print(test_dna.translate_seq())
print(test_dna.codon_usage('C'))
print(test_dna.gen_reading_frames())

for rf in test_dna.gen_reading_frames():
    print(rf)

print(test_dna.all_proteins_from_orfs())

writeTextFile('test.txt',test_dna.seq)
for rf in test_dna.gen_reading_frames():
    writeTextFile('test.txt', str(rf),'a')
