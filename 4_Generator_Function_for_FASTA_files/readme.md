1) A Generator Function that reads a Fasta file. In each iteration, the function
must return a tuple with the following format: (identifier, sequence).
Function name: FASTA_iterator(fasta_filename)

2) Given a list of FASTA files, create a function that returns a dictionary that contains
the 4 following keys with the associated values:
    - “intersection”: a set with the common identifiers found in all the files
    - “union”: a set with all the identifiers (unique) found in all the files
    - “frequency”: a dictionary with all the identifiers as keys and the number of files
    in which it appears as values (int)
    - “specific”: a dictionary with the name of the input files as keys and a set with the
    specific identifiers as values (i.e. identifiers that are exclusive in that fasta file)

Note 1: Common identifier equivalence must be case-insensitive (i.e. Code_A,code_a and
CODE_A are equivalents).

Note 2: It must use the FASTA_iterator function created in exercise 1.

Function name:
compare_fasta_file_identifiers( fasta_filenames_list )