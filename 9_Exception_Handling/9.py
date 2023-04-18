import sys
import os

def errprint(text_string):
    sys.stderr.write(text_string)
    sys.stderr.flush()

# 1. Create a new ValueError exception subclass named IncorrectSequenceLetter
class IncorrectSequenceLetter(ValueError):
    def __init__(self, invalid_chars, class_name):
        self.invalid_chars = invalid_chars
        self.class_name = class_name
    
    def __str__(self):
        # 3)The description string of the exception must be the following:
        # “The sequence item B is not found in the alphabet of class ProteinSequence”
        return f"The sequence item {', '.join(list(self.invalid_chars))} is not found in the alphabet of class {self.class_name}"


class Sequence:
    alphabet = set()

    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence.upper()

        # 4)Sequence class should raise an IncorrectSequenceLetter exception when a sequence
        # is created using an incorrect letter not found in the alphabet.
        if not set(self.__sequence).issubset(self.alphabet):
            invalid_chars = set(self.__sequence) - self.alphabet
            # 2)To create a new exception instance, it should be created with the letter not found in the
            # alphabet and the class name of the sequence. Example:
            raise IncorrectSequenceLetter(invalid_chars,self.__class__.__name__)
            
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

    # ## delete after:
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


# 5)Modify the FASTA_iterator to be able to iterate in any type of Sequence (add a new
# argument in the generator function to specify the Sequence class to use)


def FASTA_iterator(fasta_filename, SeqType):
    identifier = ''
    sequence = ''
    i = 0

    with open(fasta_filename) as file:
        for line in file:
            if '>' not in line:
                sequence += line.replace('\n','')
            else:
                if i != 0:
                    # 6) Modify the FASTA_iterator generator function to skip sequences having incorrect
                    # letters . It must capture specifically the IncorrectSequenceLetter exception. When it
                    # happens, it should print a message of error in the standard error and continue with the
                    # next sequence. Do not handle other types of exceptions, only IncorrectSequenceLetter
                    # exceptions.
                    try:
                        yield SeqType(identifier, sequence)
                    except IncorrectSequenceLetter as ISL:
                        errprint(str(ISL) + '\n')
                i += 1
                identifier = line.replace('\n','').replace('>','')
                sequence = ''
        # 6)
        try:
            yield SeqType(identifier, sequence)
        except IncorrectSequenceLetter as ISL:
            errprint(str(ISL) + '\n')


# 7) When the script is executed as a standalone application (without being imported (i.e.
# code under __main__ block), the script should read input DNA FASTA file(s) and calculate
# the length and molecular weight of their corresponding proteins (i.e. corresponding to
# ProteinSequence instances obtained after translation). The script must print the output
# to standard output or to a file. Output should be sorted by molecular weight, from
# lowest to greatest.
# python3 NIE_exercise_block2_part4.py [IN] [OUT]

def main():

    # START: ####################################################################
    # Look at argumnets 
    input_path = None
    output_path = None

    if len(sys.argv) == 1:
        # 7.1) If the script is executed without arguments, it looks for all “.fasta” or “.fa” files in the
        # current directory, and process all of them with a single sorted output. Print the results to
        # standard output.
        input_path = '.'

    elif len(sys.argv) == 2:
        # 7.2) If the script has a single argument, it corresponds to the input. If it is a directory, it
        # looks for all “.fasta” or “.fa” files in the given directory, process them and print the results
        # in standard output. If this single argument corresponds to a file (not necessarily “.fasta”, or
        # “.fa”), process it and print the results in standard output.
        input_path = sys.argv[1]

    elif len(sys.argv) == 3:
        # If the script has two arguments, the first one corresponds to the input and the
        # second one to the output file.
        input_path = sys.argv[1]
        output_path = sys.argv[2]

    else:
        print(f"Usage: \tpython3 {sys.argv[0]} [IN] [OUT] or\
        \n\tpython3 {sys.argv[0]} [IN] or\
        \n\tpython3 {sys.argv[0]}\n", end='', sep='')


    if input_path != None:  # I included this so that the script can be run 
                            # in a notebook without error, otherwise os returns
                            # an error because the path is set to None in the start

    # NEXT STEP: ####################################################################3
    # See if input is file or directory

        files = []
        if os.path.isdir(input_path):
            files = [os.path.join(input_path, f) \
                    for f in os.listdir(input_path) \
                    if f.endswith('.fasta') or f.endswith('.fa')]

        elif os.path.isfile(input_path):
            files = [input_path]

        # Go throught files

        if files:
            errprint(f'{len(files)} FASTA files found.\n')
            results = []
            seq_count = 0
            for file in files:
                for i in FASTA_iterator(file, DNASequence):
                    translated_to_protein = ProteinSequence(i.get_identifier(), i.translate())
                    ident_len_t_mw_t = [i.get_identifier(), len(translated_to_protein), translated_to_protein.get_mw()]
                    # print(i.get_identifier())
                    # print(len(translated_to_protein))
                    # print(translated_to_protein.get_mw())
                    # print(ident_len_t_mw_t)
                    results.append(ident_len_t_mw_t)
                    seq_count += 1
                errprint(f'{file} finished.\n')
            errprint(f'{seq_count} sequences found.\n')
        else:
            errprint("Couldn't load any files\n")
            sys.exit(1)

        # Sorting:
        
        errprint("Sorting the sequences...\n")
        sorted_results = sorted(results, key=lambda x: x[2], reverse=False)
        errprint("Sort process finished.\n")
        # print(sorted_results)

        # Write oupput or print

        if output_path:
            with open(output_path, 'w') as f:
                for seq in sorted_results:
                    # print(seq)
                    # for elem in seq:
                    line = seq[0] + '\t' +  str(seq[1]) + '\t' + str(seq[2]) + '\n' 
                    f.write(line)
        else:
            for seq in sorted_results:
                # print(seq)
                # for elem in seq:
                line = seq[0] + '\t' +  str(seq[1]) + '\t' + str(seq[2]) + '\n'
                print(line, sep='', end='')

        errprint("Program finished correctly.\n")
        sys.exit(0)

## Standalone definition
if __name__ == '__main__':
    main()

