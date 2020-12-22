# AALI

Amino Acid Ligand Interactions

Note: still in development but the script is usable at this point...

Python script for checking the interactions between amino acids and ligands (or amino acids and water molecules).
More specifically, it checks the distance between each amino acid in a protein model and each ligand (or each water molecule).

## Usage:

python AALI.py

Put the script AALI.py into the folder together with the input directory "PDBs" containing pdb files. Run the script with python.

The script will read PDB files from the input directory and save the results into the output directory (default "output").

For addtional options, see argument list and explanations below.

## Functionality:

Outputs a csv file (one or multiple) with all amino acids and additional information on their distances from (or interactions with) ligand(s) and/or water molecules.

Different versions of output:

- can output either the exact distance between amino acids and ligands/water or True and False values based on the given distance threshold
- it can combine all ligands into one column, reporting the smallest distance (or a single True/False value)...
- ... or it can report mulitple columns for each ligand
- in case of multiple pdb files, it can output one csv file per pdb file, or it can combine the all of the output into a single pdb file


## Options / Parameters:

output (default: exact)
  - exact (reports an exact distance)
  - tf (reports combined 1/0 values; 1 if the ligand is within the threshold distance, 0 otherwise)
  
combined (default: no)
  - yes (reports the smallest distance from all ligands in a single column, or, if output=tf, a single 1 if any of the ligands is within the threshold distance, 0 otherwise)
  - no (reports info on all ligands)
  
threshold (default: [3.5,7]) """"""""""""""""""""""TBD""""""""""""""""""""""
  - a (list of) number(s) (applicable in the case of output = tf - what will be a threshold value for which 1 will be generated, indicating a contact between amino acid and a ligand)
  
threshold_w (default: [3.5]
  - a (list of) number(s) (applicable in the case of output = tf - what will be a threshold value for which 1 will be generated, indicating a contact between amino acid and a water molecule)
  
mode (default: residue)
  - residue (checks all atoms from a given amino acid against all atoms from a given ligand and reports back the smallest one)
  - CA (checks only the distance between C alpha from amino acid against all atoms from a given ligand and reports back the smallest one)

check_ligand (default: yes)
  - no (exclude ligands)
  - yes (include columns with ligands in the output)

check_water (default: no)
  - no (exclude water)
  - yes (include columns with water molecules in the output)
  
skip (default: no)
  - no (analyze protein structure even if it has no ligand/water)
  - yes (skip protein structure if it doesn't have ligand/water)

join (default: no)
  - no (generate a csv output file for each pdb file)
  - yes (generate a single csv output file for all pdb files)
  
#### Note: expect high RAM usage if working with a lot of pdb files and selecting join='yes' (probably 1GB+ for 5000 pdb files, depends on the files themselves).

  

## General notes:

General computation time (TBD)
Testing on multiple pdb files (TBD)

Computation time CA vs all, vs with water

Water yes - only condensed mode since some protein models have a lot of water molecules

Requirements:

The script uses biopython:

https://biopython.org/

Install it with:

pip install biopython

or

conda install biopython
