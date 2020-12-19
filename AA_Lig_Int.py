# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 09:38:30 2020

@author: Muhamed Adilovic
"""

from Bio.PDB.PDBParser import PDBParser
import pandas as pd

protein = './Test/4U9S.pdb'

ID = '4U9S'

parser = PDBParser()

structure = parser.get_structure(ID, protein)

threshold = int('4')

# output = 'exact'
output = 'tf'

# combined = 'yes'
combined = 'no'

# TODO
# Add parameter skip if no H or W present

# TODO
# Add a parameter to join output of multiple proteins into a single dataframe/csv file

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
                residue_id = residue.get_id()
                hetfield = residue_id[0]
                if hetfield[0] == 'H':
                    heteros.append(residue)
    return heteros

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
                # TODO
                # Add different modes CA or residue
                distances.append(check_distance_residue(residue, H))
    return distances

def check_all_H(structure):
    """
    Parameters
    ----------
    structure : Biopython protein structure object made with PDBParser().

    Returns
    -------
    all_ligands : dataframe of all distances between each amino acid residue 
    and ligand from the protein structure, or 
    """  
    heteros = get_heteros(structure)
    # TODO
    # check if no heteros present
    #
    hetero_names = get_hetero_names(structure)
    all_ligands = pd.DataFrame()
    if output == 'exact':
        for i in range(len(heteros)):
            all_ligands[hetero_names[i]] = check_distance_protein(structure, heteros[i])
        if combined == 'yes':
            all_ligands = all_ligands.min(axis=1)
            all_ligands.name = 'distance'
    elif output == 'tf':
        for i in range(len(heteros)):
            column_id = hetero_names[i]+'_'+str(threshold)
            all_ligands[column_id] = check_distance_protein(structure, heteros[i])
        all_ligands = all_ligands < threshold
        if combined == 'yes':
            all_ligands = all_ligands.any(axis=1)
            all_ligands.name = str(threshold)
    return all_ligands

def check_all_W(structure):
    # TODO
    return

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
    AA_names, AA_nums = get_AA_names_nums(structure)
    result = pd.DataFrame()
    result['AA_name'] = AA_names
    result['AA_num'] = AA_nums
    H_distances = check_all_H(structure)
    # TODO
    # add if else for check_all_W
    return pd.concat([result,H_distances], axis=1)

# Some checks...
let_me_try = get_heteros(structure)
AA_names, AA_nums = get_AA_names_nums(structure)
hetero_names = get_hetero_names(structure)
heteros = get_heteros(structure)
H_distances = check_all_H(structure)
check_all = combine(structure)

"""
Checked against these two servers:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909059/
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3039989/
And currently, AA_Lig_Int gives the same results as these above...
"""

# output_test = open('output_test.txt', 'w')
                
                
                
                
                