# AA_Lig_Int

Note: still in development...

Python script for checking the interactions between amino acids and ligands (or amino acids and water molecules).
More specifically, it checks the distance between each amino acid in a protein model and each ligand (or each water molecule).

Planned Functionality:

Outputs a csv file with all amino acids and additional information on their interactions.

Different versions of output:

- can output either the exact distance between amino accids and ligands/water or True and False values (1 or 0) based on the given distance threshold
- it can combine all ligands into one column, reporting the smallest distance (or a single True/False value)...
- ... or it can report mulitple columns for each ligand

Planned usage:


Planned Options / Parameters:

output
  - exact (reports an exact distance)
  - tf (reports only 1/0 values)
  
combined
  - yes (reports the smallest distance from all ligands in a single column)
  - no (reports info on all ligands)
  
threshold
  - a (list of) number(s) (applicable in the case of output = tf - what will be a threshold value for which 1 will be generated, indicating a contact between amino acid and a ligand)
  
mode
  - all (checks all atoms from a given amino acid against all atoms from a given ligand and reports back the smallest one)
  - CA (checks only the distance between C alpha from amino acid against all atoms from a given ligand and reports back the smallest one)
  
water
  - no (exclude water)
  - yes (include columns with water molecules in the output)

Notes:

General computation
Computation time CA vs all
Water yes - only condensed mode since some protein models have a lot of water molecules

Requirements:

The script uses biopython:
link
Install it with:
pip install biopython
or
conda install biopython
