# AALI

Amino Acid Ligand Interactions

Python script for checking the interactions between amino acids and ligands (and/or amino acids and water molecules).
More specifically, it checks the distance between each amino acid in a protein model and each ligand (and/or each water molecule).

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
- in case of multiple pdb files, it can output one csv file per pdb file, or it can combine all of the output into a single pdb file


## Options / Parameters:

output (default: exact)
  - exact (reports an exact distance)
  - tf (reports combined True or False values; True if the ligand is within the threshold distance, False otherwise)
  
combined (default: no)
  - yes (reports the smallest distance from all ligands in a single column, or, if output=tf, a single True if any of the ligands is within the threshold distance, False otherwise)
  - no (reports info on all ligands)
  
threshold (default: [3.5,7])
  - a (list of) number(s) (applicable in the case of output = tf - what will be a threshold value for which True will be generated, indicating a contact between amino acid and a ligand)
  
threshold_w (default: [3.5])
  - a (list of) number(s) (applicable in the case of output = tf - what will be a threshold value for which True will be generated, indicating a contact between amino acid and a water molecule)
  
mode (default: residue)
  - residue (checks all atoms from a given amino acid against all atoms from a given ligand and reports back the smallest distance)
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
  
#### Note: expect high RAM usage if working with a lot of pdb files and selecting join='yes' (probably 1GB+ for ~5000 pdb files, depends on the files themselves).

  

## General notes:

General computation time is ~5 seconds per protein, although it can be significantly different, depending on the size of the protein and number of ligands, or if it has multiple models inside.

Computation time also depends on the parameters; e.g. mode:CA is faster then mode:all, although it is less precise, since it doesn't take sidechains into consideration.

Additionally, if calculating the distance from water molecules, the number of water molecules present can signifficantly affect the computation time.

## Requirements:

The script uses biopython:

https://biopython.org/

Install it with:

pip install biopython

or

conda install biopython

## Feedback

Feel free to ask for help, suggest features, or just send a feedback :)
