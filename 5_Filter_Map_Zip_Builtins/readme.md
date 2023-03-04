Note: Use the FASTA_iterator function created in last exercises.

1) Repeat the same exercises proposed in session 2 but using the FASTA_Iterator
function created in session 4 to read the FASTA files.

2) A function that, given a multiline FASTA file, returns the length of the sequence
with the maximum length

        get_max_sequence_length_from_FASTA_file( fasta_filename )

3) A function that, given a multiline FASTA file, returns the length of the sequence
with the minimum length

        get_min_sequence_length_from_FASTA_file ( fasta_filename )

4) A function that, given a FASTA file, returns a list of tuples (identifier, sequence)
corresponding to the sequence(s) with maximum length. The list must be sorted
by the identifier (case insensitive sorted).

        get_longest_sequences_from_FASTA_file( fasta_filename )

5) A function that, given a FASTA file, returns a list of tuples (identifier, sequence)
corresponding to the sequence(s) with minimum length. The list must be sorted by
the identifier (case insensitive sorted).

        get_shortest_sequences_from_FASTA_file( fasta_filename )

6) A function that, given a protein FASTA file, returns a dictionary with the molecular
weights of all the proteins in the file. The dictionary keys must be the protein
identifiers and the associated values must be a float corresponding to the molecular
weight.

        get_molecular_weights( fasta_filename )

7) A function that, given a protein FASTA file, returns a tuple with (identifier,
sequence) of the protein with the lowest molecular weight. If there are two or more
proteins having the minimum molecular weight, just return the first one.

        get_sequence_with_min_molecular_weight( fasta_filename )

8) A function that, given a protein FASTA file, returns the mean of the molecular
weights of all the proteins

        get_mean_molecular_weight( fasta_filename )