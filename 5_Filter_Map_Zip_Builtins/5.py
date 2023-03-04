# * #########################################################################################################

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
                    yield (identifier, sequence)
                i += 1
                identifier = line.replace('\n','').replace('>','')
                sequence=''
        yield (identifier, sequence)



# 1 #########################################################################################################
# 1.A ###########################

def get_proteins_ratio_by_residue_threshold(filename,residue,relative_threshold=0.03, absolute_threshold=10):
##########
# filename = 'test.fasta'
# residue = 'Y'
# relative_threshold=0.03
# absolute_threshold=10
##########
    residue = str(residue.upper())

    num_of_proteins = 0
    count_of_true_proteins = 0
    given_residue_count = 0
    total_residue_count = 1

    fd = FASTA_iterator(filename)

    for line in fd: 
        line = line[1]

        total_residue_count=0
        given_residue_count=0
        num_of_proteins += 1

        given_residue_count+=line.replace("\n","").count(residue)
        total_residue_count+=len(line.replace("\n",""))

        if (given_residue_count >= absolute_threshold) and ((given_residue_count/total_residue_count) >= relative_threshold):
            count_of_true_proteins += 1

    ratio = float(count_of_true_proteins/num_of_proteins)

    return ratio

# 1.B ###########################

def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):
##########
# filename = 'test.fasta'
# output_filename = 'test.output'
# first_n=3
# last_m=4
# #######
    output = open(output_filename, 'w')

    fd = FASTA_iterator(filename)
    for line in fd: 
        no_head = line[1]
        title = line[0]
        abv_and_count = ""
        c = 0

        first = no_head[0:first_n]
        last = no_head[-last_m:]

        for i in no_head:
            if no_head.count(i) != 0:
                abv_and_count += str(i) + ":" + str(no_head.count(i)) + ","
            no_head = no_head.replace(i,"")
        
        o_write = (title + "\t" + first + "\t" + last + "\t" + abv_and_count[:-1] + "\n")
        output.write(o_write)

    output.close()



# 2 #########################################################################################################
# A function that, given a multiline FASTA file, returns the length of the sequence
# with the maximum length
# get_max_sequence_length_from_FASTA_file( fasta_filename )

def get_max_sequence_length_from_FASTA_file( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_len = []
    for line in fd: 
        list_of_len.append(len(line[1]))
    return max(list_of_len)

# 3 #########################################################################################################
# A function that, given a multiline FASTA file, returns the length of the sequence
# with the minimum length
# get_min_sequence_length_from_FASTA_file ( fasta_filename )

def get_min_sequence_length_from_FASTA_file ( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_len = []
    for line in fd: 
        list_of_len.append(len(line[1]))
    return min(list_of_len)

# 4 #########################################################################################################
# A function that, given a FASTA file, returns a list of tuples (identifier, sequence)
# corresponding to the sequence(s) with maximum length. The list must be sorted
# by the identifier (case insensitive sorted).
# get_longest_sequences_from_FASTA_file( fasta_filename )

def get_longest_sequences_from_FASTA_file( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_tuples = []

    for line in fd:
        list_of_tuples.append(line)

    sorted_by_len_list_of_tuples = sorted(list_of_tuples, key= lambda x:len(x[1]), reverse=True)
    max = len(sorted_by_len_list_of_tuples[0][1])

    filtered_max = list(filter(lambda x: len(x[1]) == max, sorted_by_len_list_of_tuples))
    sorted_by_identifier_filtered_max = sorted(filtered_max, key=lambda x: x[0].upper())

    return sorted_by_identifier_filtered_max

# 5 #########################################################################################################
# A function that, given a FASTA file, returns a list of tuples (identifier, sequence)
# corresponding to the sequence(s) with minimum length. The list must be sorted by
# the identifier (case insensitive sorted).
# get_shortest_sequences_from_FASTA_file( fasta_filename )

def get_shortest_sequences_from_FASTA_file( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_tuples = []

    for line in fd:
        list_of_tuples.append(line)

    sorted_by_len_list_of_tuples = sorted(list_of_tuples, key= lambda x:len(x[1]), reverse=False)
    min = len(sorted_by_len_list_of_tuples[0][1])

    filtered_min = list(filter(lambda x: len(x[1]) == min, sorted_by_len_list_of_tuples))
    sorted_by_identifier_filtered_min = sorted(filtered_min, key=lambda x: x[0].upper())

    return sorted_by_identifier_filtered_min

# 6 #########################################################################################################
# A function that, given a protein FASTA file, returns a dictionary with the molecular
# weights of all the proteins in the file. The dictionary keys must be the protein
# identifiers and the associated values must be a float corresponding to the molecular
# weight.
# get_molecular_weights( fasta_filename )

def get_molecular_weights( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_tuples = []
    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}

    for line in fd:
        list_of_tuples.append(line)

    protein_weigth_dict = {}
    for tup in list_of_tuples:
        mw = 0
        for letter in tup[1]:
            mw += aminoacid_mw[letter.upper()]
        protein_weigth_dict.update({tup[0]:float(mw)})

    return protein_weigth_dict


# 7 #########################################################################################################
# A function that, given a protein FASTA file, returns a tuple with (identifier,
# sequence) of the protein with the lowest molecular weight. If there are two or more
# proteins having the minimum molecular weight, just return the first one.
# get_sequence_with_min_molecular_weight( fasta_filename )

def get_sequence_with_min_molecular_weight( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_tuples = []
    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}

    for line in fd:
        list_of_tuples.append(line)

    protein_weigth_dict = {}
    for tup in list_of_tuples:
        mw = 0
        for letter in tup[1]:
            mw += aminoacid_mw[letter.upper()]
        protein_weigth_dict.update({tup[0]:float(mw)})

    sorted_list_of_tuples = sorted(protein_weigth_dict.items(), key=lambda x:x[1], reverse=False)

    return sorted_list_of_tuples[0]


# 8 #########################################################################################################3
# A function that, given a protein FASTA file, returns the mean of the molecular
# weights of all the proteins
# get_mean_molecular_weight( fasta_filename )

def get_mean_molecular_weight( fasta_filename ):
# fasta_filename = 'test.fasta'
    fd = FASTA_iterator(fasta_filename)
    list_of_tuples = []
    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}

    for line in fd:
        list_of_tuples.append(line)

    list_with_mw = []
    for tup in list_of_tuples:
        mw = 0
        for letter in tup[1]:
            mw += aminoacid_mw[letter.upper()]
        list_with_mw.append(float(mw))

    average = sum(list_with_mw)/len(list_with_mw)

    return average

#############################################################################################################
