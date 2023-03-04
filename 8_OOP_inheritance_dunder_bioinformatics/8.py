
class Sequence:
    alphabet = set()

    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence

        if not set(self.__sequence).issubset(self.alphabet):
            invalid_chars = set(self.__sequence) - self.alphabet
            raise ValueError(f"Impossible to create instance: {invalid_chars} not possible")

    def get_identifier(self):
        return self.__identifier

    def get_sequence(self):
        return self.__sequence
        
    def get_mw(self):

        protein_weights = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        rna_weights = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
        dna_weights = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}

        if isinstance(self, ProteinSequence):
            weigths_dict = protein_weights
        elif isinstance(self, DNASequence):
            weigths_dict = dna_weights
        elif isinstance(self, RNASequence):
            weigths_dict = rna_weights
        
        return sum(weigths_dict[monomer] for monomer in self.get_sequence())

    def has_subsequence(self, sequence_obj):
        return sequence_obj.get_sequence() in self.__sequence

    # len(Sequence): should return the length of the sequence.
    def __len__(self):
        return len(self.get_sequence())

    # sequence1 == sequence2: return True if sequence strings are exactly the
    # same (without taking into account the identifiers).
    def __eq__(self, other):
        return self.get_sequence() == other.get_sequence()

    # sequence1 != sequence2: return True if sequences are different, without
    # taking into account the identifiers.
    def __ne__(self, other):
        return self.get_sequence() != other.get_sequence()

    # Sequence + Sequence: Create a new sequence object instance with their
    # sequences concatenated. Sequence object has to be of the same class as the
    # operands. It should not be applicable to different classes (i.e. ProteinSequence,
    # RNASequence). The identifiers should also be concatenated with a “+” as a glue
    # between both identifiers.
    def __add__(self, other):
        if type(self) != type(other):
            raise TypeError("Can only combine same type of sequences!")
        
        comb_ident = self.get_identifier() + "+" + other.get_identifier()
        comb_seq = self.get_sequence() + other.get_sequence()

        return self.__class__(comb_ident, comb_seq)

    # Sequence[i]: should return the sequence element at position i. Position 0
    # corresponds to the first position.
    def __getitem__(self,key):
        return self.get_sequence()[key]

    # in operator: should return a boolean if the string is a substring of the attribute
    # sequence.
    def __contains__(self, item):
        return item in self.get_sequence()

    # Comparing sequences. Implement the necessary method(s) to define how
    # sequences should be ordered. The objective is that when sorting a list of
    # sequences, they are sorted according to their molecular weight.
    def __lt__(self, other):
        return self.get_mw() < other.get_mw() 
    def __le__(self, other):
        return self.get_mw() <= other.get_mw() 
    def __gt__(self, other):
        return self.get_mw() > other.get_mw() 
    def __ge__(self, other):
        return self.get_mw() >= other.get_mw() 

    # Adapt the sequence class so that it can be used as key in a dictionary or it can be
    # added to a set. Two sequences should be considered the same object in terms of
    # set or key if they share both the identifier and the sequence.
    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))

    # ## delete after
    # def __str__(self) -> str:
    #     return self.get_identifier() + " == " + self.get_sequence()
    # ###


class ProteinSequence(Sequence):
    protein_letters = 'ACDEFGHIKLMNPQRSTVWY'
    alphabet = set(protein_letters)


class NucleotideSequence(Sequence):
    
    def translate(self):

        rna_table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
        rna_stop_codons = ['UAA', 'UAG', 'UGA']
        rna_start_codons = ['UUG', 'CUG', 'AUG']

        dna_table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}
        dna_stop_codons = ['TAA', 'TAG', 'TGA']
        dna_start_codons = ['TTG', 'CTG', 'ATG']

        if isinstance(self, DNASequence):
            codon_table_dict = dna_table
            start_list = dna_start_codons
            stop_list = dna_stop_codons
        elif isinstance(self, RNASequence):
            codon_table_dict = rna_table
            start_list = rna_start_codons
            stop_list = rna_stop_codons
        else:
            raise ValueError('Invalid sequence type')

        sequence = self.get_sequence()
        aa_list = []

        has_start_index = [sequence.find(codon) for codon in start_list \
            if codon in sequence]
        
        if has_start_index:
            start_index = min([sequence.find(codon) for codon in start_list \
                if codon in sequence])
            for i in range(start_index, len(sequence), 3):
                codon = sequence[i: i+3]
                if codon in stop_list or len(codon) != 3 :
                    break
                aa_list.append(codon_table_dict[codon])
            
            return ''.join(aa_list)

        else:
            raise ValueError('Object has no start codon')

class DNASequence(NucleotideSequence):
    dna_letters = 'GATC'
    alphabet = set(dna_letters)
    
    def transcribe(self):
        dna_to_rna = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([dna_to_rna[D] for D in self.get_sequence()])


class RNASequence(NucleotideSequence):
    rna_letters = 'GAUC'
    alphabet = set(rna_letters)

    def reverse_transcribe(self):
        rna_to_dna = {'U':'A', 'G':'C', 'C':'G', 'A':'T'}
        return ''.join([rna_to_dna[R] for R in self.get_sequence()])



# # DNA Sequences
# dna1 = DNASequence("dna1", "ATCG")
# dna2 = DNASequence("dna2", "GCTA")
# dna3 = DNASequence("dna3", "ATGCATGC")
# dna4 = DNASequence("dna4", "CGTA")
# dna5 = DNASequence("dna5", "TTGGGCGCGCGCCGCTAGCGCG")

# # RNA Sequences
# rna1 = RNASequence("rna1", "AUCG")
# rna2 = RNASequence("rna2", "CGAU")
# rna3 = RNASequence("rna3", "UAGCUAGC")
# rna4 = RNASequence("rna4", "UAUG")
# rna5 = RNASequence("rna5", "UGCGUGCGUGCGUGCGUGCGUGCGUGCGUGCGUGCG")

# # Protein Sequences
# prot1 = ProteinSequence("prot1", "MV")
# prot2 = ProteinSequence("prot2", "MMYTFGGVSQPTPY")
# prot3 = ProteinSequence("prot3", "MSKGEELFTGYVQILGHKLEYNTRDYPQIDARYYREQVKGVVVV")
# prot4 = ProteinSequence("prot4", "MKLAILFVVAVFASASA")
# prot5 = ProteinSequence("prot5", "MTRRAQMAST")


# objs_lsit = [
#     dna1, dna2, dna3, dna4, dna5, \
#     rna1, rna2, rna3, rna4, rna5, \
#     prot1, prot2, prot3, prot4, prot5]


# dna1 = DNASequence("seq1", "ATCG")
# dna_1 = DNASequence("seq_1", "ATCG")
# dna2 = DNASequence("seq2", "GCTA")
# dna3 = DNASequence("seq1", "ATCG")

# dna_dict = {dna1: "sequence 1", dna2: "sequence 2"}
# print(dna_dict[dna3])  # prints "sequence 1"

# dna_set = {dna1, dna2, dna3}
# print(len(dna_set))  # prints 2
