# 1) Define a new class named Protein, with the following definition. If
# necessary, you can define the private methods you need
# -------------------
# Protein
# -------------------
# +identifier: String
# +sequence: String
# -------------------
# +get_identifier(): string
# +get_sequence(): string
# +get_mw(): float
# +has_subsequence( Protein): boolean
# +get_length(): integer

global aminoacid_mw
aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, \
    'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19,\
    'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, \
    'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}


class Protein():

    def __init__(self, identifier, sequence):
        self.identifier = identifier
        self.sequence = sequence
    
    def get_identifier(self):
        return self.identifier

    def get_sequence(self):
        return self.sequence

    def get_mw(self):
        mw = 0
        seq = self.sequence
        for i in seq:
            mw += aminoacid_mw[i.upper()]
        return mw
    
    def has_subsequence(self, Protein):
        seq = self.sequence
        subseq = Protein.get_sequence()
        return subseq in seq

    def get_length(self):
        return (len(self.sequence))

# a = Protein(identifier="TEST", sequence="SASADASDASD")
# b = Protein("test2", "ADAS")
# c = Protein('test3', 'KTR')

# * #########################################################################################################
# 2) Modify the FASTA_iterator generator function to yield Protein objects
# instead of tuples. In each iteration, the function must yield a Protein
# Object:
# FASTA_iterator( fasta_filename )

def FASTA_iterator(fasta_filename):
    identifier = ''
    sequence = ''
    i = 0
    with open(fasta_filename) as file:
        for line in file:
            if '>' not in line:
                sequence += line.replace('\n','')
            else:
                if i != 0:
                    yield Protein(identifier, sequence)
                i += 1
                identifier = line.replace('\n','').replace('>','')
                sequence=''
        yield Protein(identifier, sequence)


# gen = FASTA_iterator('input.fasta')
# print('identifier \tlength \tmw')
# for ent in gen:
#     print(ent.get_identifier(), '', ent.get_length(), ent.get_mw(), sep= '\t')
    