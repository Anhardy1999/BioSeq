from BioStructs import DNA_Codons, RNA_Codons, NUCLEOTIDE_BASE
from collections import Counter
import random


class BioSeq:
    '''DNA sequence class. Default value: ACGT, DNA. No label
    Label: FASTA file or Database files header. Keep track of what you are working with. '''
    def __init__(self, seq = 'ACGT', seq_type = 'DNA', label = 'No Label'):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validateSeq()
        assert self.is_valid, f'Provided data does not seem to be a correct {self.seq_type} sequence.'

    # DNA toolkit methods

    def __validateSeq(self):
        ''' Ensure each randomly generated DNA sequence is a valid sequence. '''
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        ''' Returns sequence type. '''
        return self.seq_type

    def get_seq_info(self):
        ''' Returns 4 strings. Full sequence information.
        Label, Sequence, Sequence Type, and Length
        '''
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rand_seq(self, length = 10, seq_type = 'DNA'):
        ''' Generate a random sequence provided the length.
        The default sequence is the sequence that is provided during the 
        creation of this model.
        '''

        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
        for x in range(length)])
        self.__init__(seq, seq_type, 'Randomly Generated Sequence')

    def nucleotide_freq(self):
        '''Count the nucleotides in a given sequence. Returns a dictionary. '''
        return dict(Counter(self.seq))

    def transcription(self):
        ''' DNA -> RNA transcription. Replacing Thymine with Uracil. '''
        
        if self.seq_type == 'DNA':
            return self.seq.replace('T', 'U')

        return 'Not a DNA sequence'

    def reverse_complement(self):
        ''' Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string
        '''

        if self.seq_type == 'DNA':
            mapping = str.maketrans('ACGT', 'TGCA')
        
        else:
            mapping = str.maketrans('ACGU', 'UGCA')
        
        return self.seq.translate(mapping)[::-1]    

    def gc_content(self):
        ''' GC Content in DNA/RNA sequence.
        Important when  designing primers or cloning projects. High GC contents can result
        in stronger DNA binding with the primer or the other complementary strand.  Lower GC-content
        (>40%) is ideal. Also determines the DNA melting termperature to separate the single strands for efficient
        PCR reactions 
        '''
        return round(((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100), 5)

    def gc_content_subset(self, k = 20):
        ''' GC Content in a DNA/RNA sub-sequence length k. k = 20 by default '''
        res = []
        for i in range(0, len(self.seq) - k + 1, k): # looping through the sequence at a specific window size.
            subseq = self.seq [i:i + k]
            res.append(round(((subseq.count('C') + subseq.count('G')) / len(subseq) * 100), 5))
        
        return res

    def translate_seq(self, init_pos = 0):
        ''' Translates the sequence into an amino acid sequence 

        init_pos: For adjusting the reading frames. By default it is 0. 
        '''
        if self.seq_type == 'DNA':
            return ''.join([DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2,3)])
        
        elif self.seq_type == 'RNA':
            return ''.join([RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2,3)])

    def codon_usage(self, aminoacid):
        ''' Calculates the frequency of each codon encoding for a given amino acid in a sequence '''
        tmpList = []
        if self.seq_type == 'DNA':
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        elif self.seq_type == 'RNA':
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWeight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWeight, 2)

        return freqDict

    def gen_reading_frames(self):
        ''' Generate the six reading frames of sequence, including the reversed complement '''
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        
        temp_seq = BioSeq(self.reverse_complement(), self.seq_type)
        frames.append(temp_seq.translate_seq(0))
        frames.append(temp_seq.translate_seq(1))
        frames.append(temp_seq.translate_seq(2))

        del temp_seq

        return frames

    def proteins_from_rf(self, aa_seq):
        ''' Compute all possible proteins in an AA seq and return a list of possible proteins
        Takes a list of Proteins'''
        current_prot = [] # Will contain the current proteins
        proteins = [] # contains all of proteins found in AA sequence
        for aa in aa_seq:
            if aa == '_':
                # STOP accumulating Amino Acids if _. Stop was found
                if current_prot: # checks if length is more than 1 basically
                    for p in current_prot: # make sure that if we have multiple stop proteins we can keep going
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating AA if M - START was found
                if aa == 'M':
                    current_prot.append('')
                for i in range(len(current_prot)):
                    current_prot[i] += aa

        return proteins

    def all_proteins_from_orfs(self, startReadPos = 0, endReadPos = 0, ordered = False):
        ''' Compute all possible proteins for all open reading frames '''

        if endReadPos > startReadPos:
            temp_seq = BioSeq(self.seq[startReadPos:endReadPos])
            rfs = temp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()
        
        res = []
        for rf in rfs: 
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        
        if ordered:
            return sorted(res, key=len, reverse = True)
        return res

    def hamming_distance(self, str1, str2): 
        ''' Returns the hamming distance of two provided sequences '''
        if len(str1) == len(str2):
            return sum([1 for i in range(len(str1)) if str1[i]!=str2[i]])     
