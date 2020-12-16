# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 09:38:30 2020

@author: Muhamed Adilovic
"""

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection

def check_distance(a, b):
    return a-b

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

let_me_try = get_heteros(structure)