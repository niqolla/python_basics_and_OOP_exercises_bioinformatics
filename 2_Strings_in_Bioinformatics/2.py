# 1) Given a multi-line protein FASTA file (stored in a file with path defined filename), returns a
# float corresponding to the ratio of proteins in the fasta file having a relative frequency higher
# or equal than a given threshold provided as an argument named “relative_threshold” and
# having an absolute frequency of the same residue higher or equal than a given threshold
# provided as an argument named “absolute_threshold” for a given residue. The function
# should be named as follows, with the same arguments definition:

def get_proteins_ratio_by_residue_threshold(filename,residue,relative_threshold=0.03, absolute_threshold=10):

    residue = str(residue.upper())

    num_of_proteins = 0
    count_of_true_proteins = 0
    given_residue_count = 0
    total_residue_count = 1
    # I'm setting total_residue_count to 1 before the start of the function so that 
    # in the first itteration there won't be an issue with diving by 0 (if absolute_threshold = relative_threshold = 0)
    # after the if, the total_residue_count is set 0

    with open(filename, "r") as fd:
        for line in fd: 
            if ">" in line:

                if (given_residue_count >= absolute_threshold) and ((given_residue_count/total_residue_count) >= relative_threshold):
                    count_of_true_proteins += 1

                total_residue_count=0
                given_residue_count=0
                num_of_proteins += 1

            else:
                
                given_residue_count+=line.replace("\n","").count(residue)
                total_residue_count+=len(line.replace("\n",""))

        if (given_residue_count >= absolute_threshold) and ((given_residue_count/total_residue_count) >= relative_threshold):
            count_of_true_proteins += 1
                
    ratio = float(count_of_true_proteins/num_of_proteins)

    return ratio


# 2) Given a protein FASTA file (filename), save on a output file named output_filename the
# protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency
# in the protein of all the aminoacids found in the protein (the aminoacids that do not appear
# in the protein should not be shown). The fields must be separated by a tabulator, and one
# protein by line.  

def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):

    no_head=""
    title = ""
    abv_and_count = ""
    c = 0

    output = open(output_filename, 'w')

    with open(filename, "r") as fd:
        for line in fd: 
            if ">" in line:
                if c != 0:
                    output.write( title + "\t" + no_head[0:first_n] + "\t" + no_head[-last_m:] + "\t")
                    for i in no_head:
                        if no_head.count(i) != 0:
                            abv_and_count += str(i) + ":" + str(no_head.count(i)) + ","
                        no_head = no_head.replace(i,"")
                    output.write(abv_and_count[:-1] +"\n")
                
                no_head=""
                abv_and_count = ""
                c+=1
                title = line.replace(">","").replace("\n","")

            else:
                no_head+=line.replace("\n","")
        output.write(title + "\t" + no_head[0:first_n] + "\t" + no_head[-last_m:] + "\t")
        for i in no_head:
            if no_head.count(i) != 0:
                abv_and_count += str(i) + ":" + str(no_head.count(i)) + ","
            no_head = no_head.replace(i,"")
        output.write(abv_and_count[:-1] + "\n")
    output.close()
