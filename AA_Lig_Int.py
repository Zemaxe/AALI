# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 09:38:30 2020

@author: Muhamed Adilovic
"""

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection

protein = './Test/4U9S.pdb'

ID = '4U9S'

parser = PDBParser()

structure = parser.get_structure(ID, protein)

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

let_me_try = get_heteros(structure)

for H in let_me_try:
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_nonAA(residue):
                    continue
                else:
                    print(residue.get_resname())
                    print(check_distance_CA(residue, H))
                
                
                
                
                
                
                
                