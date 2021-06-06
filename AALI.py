# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 09:38:30 2020

@author: Muhamed Adilovic
"""

import os
import argparse
import pandas as pd
from Bio.PDB.PDBParser import PDBParser

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", choices=["exact", "tf"], default="exact",
                        help = "set output type: exact distance or true/false")
    parser.add_argument("-in", "--in_folder", default="PDBs",
                        help = "set the folder containing pdb files")
    parser.add_argument("-t", "--threshold", type=int, default=4,
                        help = "set the threshold value for generating True/False output")
    parser.add_argument("-out", "--out_folder", default="output",
                        help = "set the folder which will contain the results")
    parser.add_argument("-c", "--combined", choices=["yes", "no"], default="no",
                        help = "choose whether to show all distances or the smallest ones")
    parser.add_argument("-m", "--mode", choices=["CA", "residue"], default="residue",
                        help = "choose whether to check distance against CA or \
                            all atoms in amino acid and pick the smallest one")
    parser.add_argument("-chl", "--check_ligand", choices=["yes", "no"], default="yes",
                        help = "choose whether to check amino acids against ligands")
    parser.add_argument("-chw", "--check_water", choices=["yes", "no"], default="no",
                        help = "choose whether to check amino acids against water")
    parser.add_argument("-s", "--skip", choices=["yes", "no"], default="no",
                        help = "choose whether to analyze the pdb file even if \
                            it doesn't have water or ligands")
    parser.add_argument("-j", "--join", choices=["yes", "no"], default="no",
                        help = "choose whether to join all output csv files into one")
    args = parser.parse_args()
    return args

threshold_w = int('4')

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
    AA_model = []
    AA_chain = []
    for model in structure:
        for chain in model:
            for residue in chain.get_list():
                if is_nonAA(residue):
                    continue
                AA_names.append(residue.get_resname())
                AA_nums.append(residue.get_id()[1])
                AA_model.append(model._id)
                AA_chain.append(chain._id)               
    return AA_names, AA_nums, AA_model, AA_chain
                
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
                if args.mode == 'residue':
                    distances.append(check_distance_residue(residue, H))
                elif args.mode == 'CA':
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
    hetero_names = get_hetero_names(structure)
    all_ligands = pd.DataFrame()
    # handling a case if output = 'exact'
    if args.output == 'exact':
        for i in range(len(heteros)):
            all_ligands[hetero_names[i]] = check_distance_protein(structure, heteros[i])
        if args.combined == 'yes':
            all_ligands = all_ligands.min(axis=1)
            all_ligands.name = 'ligand_distance'
    # handling a case if output = 'tf'
    elif args.output == 'tf':
        for i in range(len(heteros)):
            column_id = hetero_names[i]+'_'+str(args.threshold)
            all_ligands[column_id] = check_distance_protein(structure, heteros[i])
        all_ligands = all_ligands < args.threshold
        if args.combined == 'yes':
            all_ligands = all_ligands.any(axis=1)
            all_ligands.name = 'ligand_'+str(args.threshold)
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
    all_water = pd.DataFrame()
    # handling a case if output = 'exact'
    if args.output == 'exact':
        for i in range(len(water)):
            all_water[str(i)] = check_distance_protein(structure, water[i])
        all_water = all_water.min(axis=1)
        all_water.name = 'water_distance'
    # handling a case if output = 'tf'
    elif args.output == 'tf':
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
    None : in cases where the protein is supposed to be skipped
    """
    # TODO
    # Add support to generate an empty column if no ligand/water is present
    # and skip == no
    
    AA_names, AA_nums, AA_model, AA_chain = get_AA_names_nums(structure)
    result = pd.DataFrame()
    result['AA_name'] = AA_names
    result['AA_num'] = AA_nums
    result['model'] = AA_model
    result['chain'] = AA_chain
    heteros = get_heteros(structure)
    water = get_water(structure)
    # handle different types of checks...
    if args.check_ligand == 'yes' and args.check_water == 'yes':
        if args.skip == 'yes' and len(heteros) == 0 and len(water) == 0:
            return
        H_distances = check_all_H(structure)
        W_distances = check_all_W(structure)
        return pd.concat([result,H_distances,W_distances], axis=1)
    elif args.check_ligand == 'yes':
        if args.skip == 'yes' and len(heteros) == 0:
            return
        H_distances = check_all_H(structure)
        return pd.concat([result,H_distances], axis=1)
    elif args.check_water == 'yes':
        if args.skip == 'yes' and len(water) == 0:
            return
        W_distances = check_all_W(structure)
        return pd.concat([result,W_distances], axis=1)
    else:
        print('Please use "yes" for at least one of the two: chek_water or check_ligand.')
        return

def main(args):
    if not os.path.exists(args.out_folder):
        os.mkdir(args.out_folder)
    
    proteins = []
    
    for file in os.listdir(args.in_folder):
        if file.endswith('.pdb'):
            proteins.append(file)
    
    if args.join == 'yes':
        result_joined = pd.DataFrame()
    
    for protein in proteins:
        ID = protein.replace('.pdb', '')
        parser = PDBParser()
        protein_path = args.in_folder+'/'+protein
        structure = parser.get_structure(ID, protein_path)
        
        result = combine(structure)
        if result is None:
            print('No ligands and/or water present in ', ID)
            continue
        if args.join == 'no':
            result_path = args.out_folder+'/'+ID+'.csv'
            result.to_csv(result_path)
        elif args.join == 'yes':
            result_joined = pd.concat([result_joined,result])
    
    if args.join == 'yes':
        result_joined_path = args.out_folder+'/'+'AALI_contacts.csv'
        result_joined.to_csv(result_joined_path)

if __name__ == '__main__':
    args = parse_args()
    main(args)

                
                
                
                
                