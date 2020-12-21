# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 09:38:30 2020

@author: Muhamed Adilovic
"""

from Bio.PDB.PDBParser import PDBParser
import pandas as pd
import os

threshold = int('4')

threshold_w = int('4')

mode = 'residue'
# mode = 'CA'

output = 'exact'
# output = 'tf'

# combined = 'yes'
combined = 'no'

skip = 'yes'
# skip = 'no'

check_water = 'yes'
# check_water = 'no'

check_ligand = 'yes'
# check_ligand = 'no'

in_folder = 'PDBs'
out_folder = 'output'

join = 'no'
# join = 'yes'

# TODO
# Implement handling multiple chains
# E.g. if multiple chains add chain column in the output

# TODO
# Implement support for verbose
# Printing out the percentage of work done

# TODO
# Implement adding a column with protein name in case of join == 'yes'

def get_heteros(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    heteros : a list of hetero residues in a protein model.
    """
    heteros = []

    for model in structure:
        for chain in model:
            for residue in chain.get_list():
                if is_hetero(residue):
                    heteros.append(residue)
    return heteros

def get_water(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    water : a list of water molecules in a protein model.
    """
    water = []

    for model in structure:
        for chain in model:
            for residue in chain.get_list():
                if is_water(residue):
                    water.append(residue)
    return water

def get_hetero_names(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    hetero_names : a list of hetero residue names in a protein model.
    """
    hetero_names = []
    for model in structure:
        for chain in model:
            for residue in chain.get_list():
                if is_hetero(residue):
                    hetero_names.append(residue.get_resname())
    return hetero_names

def get_AA_names_nums(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    AA_names : a list of amino acid names in a protein model.
    AA_nums : a list of animo acid number IDs in a protein model.
    """
    AA_names = []
    AA_nums = []
    for model in structure:
        for chain in model:
            for residue in chain.get_list():
                if is_nonAA(residue):
                    continue
                AA_names.append(residue.get_resname())
                AA_nums.append(residue.get_id()[1])
    return AA_names, AA_nums
                
def is_nonAA(residue):
    """
    Parameters
    ----------
    residue : a residue from a protein structure object made with PDBParser().

    Returns
    -------
    Boolean
        True if residue is hetero or water, False otherwise.
    """
    residue_id = residue.get_id()
    hetfield = residue_id[0]
    return (hetfield[0] == 'H') or (hetfield[0] == 'W')

def is_hetero(residue):
    """
    Parameters
    ----------
    residue : a residue from a protein structure object made with PDBParser().

    Returns
    -------
    Boolean
        True if residue is hetero, False otherwise.
    """
    residue_id = residue.get_id()
    hetfield = residue_id[0]
    return hetfield[0] == 'H'

def is_water(residue):
    """
    Parameters
    ----------
    residue : a residue from a protein structure object made with PDBParser().

    Returns
    -------
    Boolean
        True if residue is water, False otherwise.
    """
    residue_id = residue.get_id()
    hetfield = residue_id[0]
    return hetfield[0] == 'W'

def check_distance_CA(AA, H):
    """
    Parameters
    ----------
    AA : amino acid from a protein structure object made with PDBParser().
    H  : hetero residue from a protein structure object made with PDBParser().

    Returns
    -------
    distance : the distance between CA of AA and H's closest atom.
    """    
    distance = []
    for atom in H:
        distance.append(AA['CA']-atom)
    return min(distance)

def check_distance_residue(AA, H):
    """
    Parameters
    ----------
    AA : amino acid from a protein structure object made with PDBParser().
    H  : hetero residue from a protein structure object made with PDBParser().

    Returns
    -------
    distance : the smallest distance between the two residues
    (includes all atoms in calculation).
    """    
    distance = []
    for atom_AA in AA:
        for atom_H in H:
            distance.append(atom_AA-atom_H)
    return min(distance)

def check_distance_protein(structure, H):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().
    H  : hetero residue from a protein structure object made with PDBParser().

    Returns
    -------
    distances : the list of smallest distances between each amino acid from 
    protein structure and a given ligand H.
    """  
    distances = []
    for model in structure:
        for chain in model:
            for residue in chain.get_list():
                if is_nonAA(residue):
                    continue
                if mode == 'residue':
                    distances.append(check_distance_residue(residue, H))
                elif mode == 'CA':
                    distances.append(check_distance_CA(residue, H))                
    return distances

def check_all_H(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    all_ligands : dataframe of all distances between each amino acid residue 
    and ligand from the protein structure.
    """  
    heteros = get_heteros(structure)
    # handling a case with no ligands present
    if len(heteros) == 0:
        print(ID, ': no ligands present.')
        if skip == 'yes':
            return
    hetero_names = get_hetero_names(structure)
    all_ligands = pd.DataFrame()
    # handling a case if output = 'exact'
    if output == 'exact':
        for i in range(len(heteros)):
            all_ligands[hetero_names[i]] = check_distance_protein(structure, heteros[i])
        if combined == 'yes':
            all_ligands = all_ligands.min(axis=1)
            all_ligands.name = 'ligand_distance'
    # handling a case if output = 'tf'
    elif output == 'tf':
        for i in range(len(heteros)):
            column_id = hetero_names[i]+'_'+str(threshold)
            all_ligands[column_id] = check_distance_protein(structure, heteros[i])
        all_ligands = all_ligands < threshold
        if combined == 'yes':
            all_ligands = all_ligands.any(axis=1)
            all_ligands.name = 'ligand_'+str(threshold)
    return all_ligands

def check_all_W(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    all_water : dataframe of all distances between each amino acid residue 
    and water molecules from the protein structure.
    """  
    water = get_water(structure)
    # handling a case with no water present
    if len(water) == 0:
        print(ID, ': no water molecules present.')
        if skip == 'yes':
            return
    all_water = pd.DataFrame()
    # handling a case if output = 'exact'
    if output == 'exact':
        for i in range(len(water)):
            all_water[str(i)] = check_distance_protein(structure, water[i])
        all_water = all_water.min(axis=1)
        all_water.name = 'water_distance'
    # handling a case if output = 'tf'
    elif output == 'tf':
        for i in range(len(water)):
            all_water[str(i)] = check_distance_protein(structure, water[i])
        all_water = all_water < threshold_w
        all_water = all_water.any(axis=1)
        all_water.name = 'water_'+str(threshold_w)
    return all_water

def combine(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    dataframe : containing every amino acid name and ID number, toegether with
    all distances between each amino acid residue and ligand from the protein 
    structure.
    """
    # TODO
    # Add support for skip entire protein if both no ligand and no water molecules are present
    # and skip == yes

    # TODO
    # Add support to generate an empty column if no ligand/water is present
    # and skip == no
    
    AA_names, AA_nums = get_AA_names_nums(structure)
    result = pd.DataFrame()
    result['AA_name'] = AA_names
    result['AA_num'] = AA_nums
    # handle different types of checks...
    if check_ligand == 'yes' and check_water == 'yes':
        H_distances = check_all_H(structure)
        W_distances = check_all_W(structure)
        return pd.concat([result,H_distances,W_distances], axis=1)
    elif check_ligand == 'yes':
        H_distances = check_all_H(structure)
        return pd.concat([result,H_distances], axis=1)
    elif check_water == 'yes':
        W_distances = check_all_W(structure)
        return pd.concat([result,W_distances], axis=1)
    else:
        print('Please use "yes" for at least one of the two: chek_water or check_ligand.')
        return
    
if not os.path.exists(out_folder):
    os.mkdir(out_folder)

proteins = []

for file in os.listdir(in_folder):
    if file.endswith('.pdb'):
        proteins.append(file)

if join == 'yes':
    result_joined = pd.DataFrame()

for protein in proteins:
    ID = protein.replace('.pdb', '')
    parser = PDBParser()
    protein_path = in_folder+'/'+protein
    structure = parser.get_structure(ID, protein_path)
    
    result = combine(structure)
    if join == 'no':
        result_path = out_folder+'/'+ID+'.csv'
        result.to_csv(result_path)
    elif join == 'yes':
        result_joined = pd.concat([result_joined,result])

if join == 'yes':
    result_joined_path = out_folder+'/'+'AALI_contacts.csv'
    result_joined.to_csv(result_joined_path)

# # Some checks...
# let_me_try = get_heteros(structure)
# AA_names, AA_nums = get_AA_names_nums(structure)
# hetero_names = get_hetero_names(structure)
# heteros = get_heteros(structure)
# water = get_water(structure)
# H_distances = check_all_H(structure)
# check_all = combine(structure)
# check_water = check_all_W(structure)



"""
Checked against these two servers:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909059/
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3039989/
And currently, AA_Lig_Int gives the same results as these above...
"""

# output_test = open('output_test.txt', 'w')
                
                
                
                
                