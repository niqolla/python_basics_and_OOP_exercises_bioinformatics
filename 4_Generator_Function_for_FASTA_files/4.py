# 1 #############################################################################

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


# 2 #############################################################################

def compare_fasta_file_identifiers(fasta_filenames_list):
    GRAND_DICT = {}

#### create a dictionary from each file:
#### { filename : indetifiers }

    for file in fasta_filenames_list:
        filename = file

        # note 2: it must use the FASTA_iterator function created in exercise 1
        file = FASTA_iterator(file)
        
        identifiers_in_file = []
        for i in file:
            # note 1: common identifier equivalence must be case-insensitive
            # here all identifiers are passed as upper
            identifiers_in_file.append(str(i[0]).upper())
    
        GRAND_DICT.update({str(filename) : identifiers_in_file})

######## creating a list of all sets

        sets = []
        for i in GRAND_DICT:
            sets.append(set(GRAND_DICT[i]))

######## “intersection”: a set with the common identifiers found in all the files

        inter_v = sets[0]
        for i in sets:
            inter_v = inter_v.intersection(i)

######## “union”: a set with all the identifiers (unique) found in all the files

        union_v = sets[0]
        for i in sets:
            union_v = union_v.union(i)

######## “frequency”: a dictionary with all the identifiers as keys 
######## and the number of files in which it appears as values (int)

        freq_v = {}
        unique_v = []
        for identifier in union_v:
            count = 0
            for i in GRAND_DICT:
                if identifier in GRAND_DICT[i]:
                    count += 1 
            if count==1:
                unique_v.append(identifier)
            freq_v.update({identifier: count})

######## “specific”: a dictionary with the name of the input files as keys 
######## and a set with the specific identifiers as values (i.e. 
######## identifiers that are exclusive in that fasta file)

        unique_file_identifiers = {}

        for n in GRAND_DICT:
            set_of_unique = set()
            
            for i in unique_v:    
                if i in GRAND_DICT[n]:
                    set_of_unique.add(i)

            if len(set_of_unique) > 0:
                unique_file_identifiers.update({n:set_of_unique})

        final_dict = {'intersection': inter_v, 'union': union_v, 'frequency': freq_v, 'specific':unique_file_identifiers}

    return final_dict
