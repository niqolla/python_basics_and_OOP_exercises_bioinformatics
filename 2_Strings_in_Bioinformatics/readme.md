Input is a multiline FASTA file.

1) Given a multi-line protein FASTA file (stored in a file with path defined filename), returns a
float corresponding to the ratio of proteins in the fasta file having a relative frequency higher
or equal than a given threshold provided as an argument named “relative_threshold” and
having an absolute frequency of the same residue higher or equal than a given threshold
provided as an argument named “absolute_threshold” for a given residue. The function
should be named as follows, with the same arguments definition:

    get_proteins_ratio_by_residue_threshold(filename,
        residue,relative_threshold=0.03, absolute_threshold=10)


2) Given a protein FASTA file (filename), save on a output file named output_filename the
protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency
in the protein of all the aminoacids found in the protein (the aminoacids that do not appear
in the protein should not be shown). The fields must be separated by a tabulator, and one
protein by line.

    print_sequence_summary(filename,
        output_filename,
        first_n=10,
        last_m=10)

Example:
    Input: 
    >PROT1
    EFTRPTSTWSAAALMTRSSSTRWSPD
    >PROT2
    SSTPLRRSTPAWEEFGLMCCDPRS
    >PROT3
    ATRSLEWKSTPW

    Output:
    PROT1 EFT RWSPD E:1,F:1,T:5,R:3,P:2,S:6,W:2,A:3,L:1,M:1,D:1
    PROT2 SST CDPRS S:4,T:2,P:3,L:2,R:3,A:1,W:1,E:2,F:1,G:1,M:1,C:2,D:1
    PROT3 ATR KSTPW A:1,T:2,R:1,S:2,L:1,E:1,W:2,K:1,P:1