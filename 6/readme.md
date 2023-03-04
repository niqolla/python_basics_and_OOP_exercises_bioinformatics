Create a python function that calculates the mean of the minimum distance between any two residues
pairs found in the same chain of a PDB. The script, when executed by command line, should output in
standard output the mean distance for each chain (with 4 decimal positions). The python script should use
a single argument corresponding to the PDB file path to use. This command line argument is optional. If
the PDB file path is not defined, read the PDB file from standard input.
The function must return a dictionary with chains as keys and mean minimum distances as values. It uses
a single argument which specifies the path of the PDB file If the argument pdb_file_path is None, read the
PDB file from standard input.


    calculate_pdb_chain_mean_minimum_distances(pdb_file_path)


When the file is imported as a module, it should not execute the function. The function should only be
called when the script is executed by command line.


Sample output in command line:
A: 22.7400
B: 20.4224
N: 23.9393
F: 23.9730
J: 23.4187