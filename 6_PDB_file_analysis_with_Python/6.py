# * #########################################################################################################
# Create a python function that calculates the mean of the minimum distance between any two residues
# pairs found in the same chain of a PDB. The script, when executed by command line, should output in
# standard output the mean distance for each chain (with 4 decimal positions). The python script should use
# a single argument corresponding to the PDB file path to use. This command line argument is optional. If
# the PDB file path is not defined, read the PDB file from standard input.
# 
# The function must return a dictionary with chains as keys and mean minimum distances as values. It uses
# a single argument which specifies the path of the PDB file If the argument pdb_file_path is None, read the
# PDB file from standard input.
# 
# calculate_pdb_chain_mean_minimum_distances(pdb_file_path)
# 
# When the file is imported as a module, it should not execute the function. The function should only be
# called when the script is executed by command line.
# * #########################################################################################################

import sys

def euclidian_distance_3d(p,q):
    d = ( (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 )**(1/2)
    return d


# pdb_file_path = '1gcn.pdb'
# pdb_file_path = '1a3n.pdb'
# pdb_file_path = 'test.pdb'
# pdb_file_path = 'simplest_test.pdb'
# pdb_file_path = 'test2.pdb'


def calculate_pdb_chain_mean_minimum_distances(pdb_file_path=None):
    # CREATING A LIST OF ALL THE ENTERIES IN THE PDB 
    # WITH 'ATOM' IN THE FIRST FIELD

    ATOMS_GRAND_LIST = []
    
    if pdb_file_path is None:
        file = sys.stdin
        for line in file:
            line_list = line.split(' ')
            while '' in line_list:
                line_list.remove('')
            while '\n' in line_list:
                line_list.remove('\n')
            if 'ATOM' in line_list[0]:
                ATOMS_GRAND_LIST.append(line_list)
    else:
        with open(pdb_file_path, 'r') as file:
            for line in file:
                line_list = line.split(' ')
                while '' in line_list:
                    line_list.remove('')
                while '\n' in line_list:
                    line_list.remove('\n')
                if 'ATOM' in line_list[0]:
                    ATOMS_GRAND_LIST.append(line_list)


    # GETTING A SET OF THE CHAINS IN THE PDB FILE
    chains = set()
    for elem in ATOMS_GRAND_LIST:
        chains.add(elem[4])
    chains = list(chains)

    # CREATING A DICTONARY WITH STRUCUTRE:
    #   {
    #   CHAIN  
    #       :
    #       [ LIST (OF RESSIDUES) OF
    #           [ LISTS (OF ATOMS)
    #               (OF COORNIATES FOR EACH ATOMS AS TUPLE (X,Y,Z) ) 
    #           ]
    #       ] 
    #   }
    #
    dict_chian_coor = {}
    ressidue = 0
    ressidue_coor = []

    for chain in chains:
        ressidue = 0
        ressidue_coor = []

        for elem in ATOMS_GRAND_LIST:

            if elem[4] == chain:

                if elem[2] == 'N':
                    if ressidue > 0 : ressidue_coor.append(coordinates_in_residue)
                    ressidue += 1
                    coordinates_in_residue = []

                coordinates = (float(elem[6]), float(elem[7]), float(elem[8]))
                coordinates_in_residue.append(coordinates)

        if ressidue > 0 : ressidue_coor.append(coordinates_in_residue)

        dict_chian_coor.update({chain:ressidue_coor})


    return_dict = {}

    for chain in dict_chian_coor:
        # for every chain
        # print(chain)
        general_dict = {}
        list_of_eu_D = []
        for residue in range(len(dict_chian_coor[chain])):
            # for all residues in a chain
            # print('  ',residue)
            residue_residue_minimum = {}

            for atom1 in range(len(dict_chian_coor[chain][residue])):
                # for atom in residue

                # print('  ',dict_chian_coor[chain][residue][atom1])
                # atom_one_residue_minimum = {}

                # loopin through all the other residues
                for residue2 in range(len(dict_chian_coor[chain])):
                    if residue2 == residue: continue
                        # go throught all the other residues of that same chain
                        # if it's the same residue, don't calulate

                    # print('     Residue: ',residue2, 'ATOMS: ', end='')
                    list_of_d_for_one_atom_and_one_residue = []

                    # looping throught all atoms of residue2
                    for atom2 in range(len(dict_chian_coor[chain][residue2])):
                        # calculate the distance between that one atom (atom1) 
                        # and all the atoms of one other residue (residue2)
                        d = euclidian_distance_3d(dict_chian_coor[chain][residue][atom1],\
                            dict_chian_coor[chain][residue2][atom2])

                        list_of_d_for_one_atom_and_one_residue.append(d)
                        # add results to a list: atom1  <----> all atoms of another residue
                        # print(' ', atom2, d, end='')

                    # print()
                    # atom_one_residue_minimum.update({residue2 : min(list_of_d_for_one_atom_and_one_residue)})

                    if residue2 not in residue_residue_minimum or min(list_of_d_for_one_atom_and_one_residue) < residue_residue_minimum[residue2]:
                        residue_residue_minimum.update({residue2 : min(list_of_d_for_one_atom_and_one_residue)})
                        # min(list_of_d_for_one_atom_and_one_residue) = for a given atom1, the distance 
                        # with the closest atom2 from another residue (residue2) in the first itteration 
                        # the atom1 the distance for atom1 to all the other residues is set,
                        # in the further itterations, if the other atoms1 are closer to the given residues 
                        # then the inital atom1, then the value is updated

            # general_dict.update({residue:residue_residue_minimum})
            for i in residue_residue_minimum:
                list_of_eu_D.append(residue_residue_minimum[i])
                # all the minimum distance values are placed into a list

        average = sum(list_of_eu_D)/len(list_of_eu_D)
        # the average is calculated

        return_dict.update({chain:average})
        # dictonary is created -> chain : average

    for chain in return_dict:
        print(chain, ': ', return_dict[chain] ,sep='')
    
    return return_dict

if __name__ == '__main__':
    if len(sys.argv) == 1:
        calculate_pdb_chain_mean_minimum_distances()
    elif len(sys.argv) == 2:
        file = sys.argv[1]
        calculate_pdb_chain_mean_minimum_distances(file)
    else:
        print("ERROR - Usage:\tpython 6.py [pdb_file_path.pdb]\n\
        OR:\tcat [pdb_file_path.pdb] | python 6.py")
