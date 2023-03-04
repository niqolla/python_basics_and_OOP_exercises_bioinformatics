################################################################################################################################################

def calculate_aminoacid_frequencies(fasta_filename, subsequences_filename, number_of_repetitions, output_filename):

    ## creating dictionary of {>title:sequence} and counting number of enteties
    seq = ""
    title=""
    name_seq_dict = {}
    input = fasta_filename
    with open(input) as file:
        for line in file:
            if ">" in line:
                name_seq_dict.update({title:seq})
                seq = ""
                title=line.replace("\n","").replace('>','')
            else:
                seq += line.replace("\n","")
        name_seq_dict.update({title:seq})
    
    del name_seq_dict['']
    total_num_of_proteins = len(name_seq_dict)

    ## creating a list of fragments and counting them
    fragments_file = subsequences_filename
    list_of_fragments = []
    with open(fragments_file) as file:
        for line in file:
            list_of_fragments.append(line.replace("\n",""))
    while "" in list_of_fragments:
        list_of_fragments.remove("")
    # I included this while loop because without the script stores '' for empty lines, and counts them too
    total_num_of_fragments = len(list_of_fragments)

    # loops -> fragments -> proteins
    threshold = number_of_repetitions
    dict_fragment_number_of_proteins={}
    for n in list_of_fragments:
    # looping throught the fragments
        number_of_proteins_that_qualify = 0
        for i in name_seq_dict:
        # for each fragment looping through every entry in the fasta file
            if name_seq_dict[i].count(n) >= threshold:
            # if number of counts of fragment for a given protein is equal or higher than the treshold
            # add one more protein that qualifes for passing the treshold
                number_of_proteins_that_qualify += 1
        if number_of_proteins_that_qualify > 0:
        # if any proteins have > 0 entries for that fragment, store in dict
            dict_fragment_number_of_proteins.update({n:number_of_proteins_that_qualify})

    sorted_list_from_dict_fragment_number_of_proteins = sorted(dict_fragment_number_of_proteins.items(), key=lambda x:x[1], reverse=1)
    # creating a list of tuples = dict_fragment_number_of_proteins.items()
    # then sorting for any tuple x by the value x[1] ,ei. by the values (not the keys)

    filename = output_filename
    output = open(filename, 'w')
    output.write("#Number of proteins:" + f"{total_num_of_proteins:>21}" + '\n')
    output.write("#Number of subsequences:" + f"{total_num_of_fragments:>17}" + '\n')
    output.write("#Subsequence proportions:" + '\n')

    for tuple_x in sorted_list_from_dict_fragment_number_of_proteins:
        dist = 21-len(tuple_x[0])
        output.write(tuple_x[0] + f"{tuple_x[1]:>{dist}}"+ f"{tuple_x[1]/total_num_of_proteins:>20.4f}" + '\n')
    output.close()

################################################################################################################################################

