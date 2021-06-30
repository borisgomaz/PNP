import os
import Bio
from Bio import PDB
io = Bio.PDB.PDBIO()
for filename in os.listdir('.'):
    if filename.endswith('.pdb'):
    
        pdb = Bio.PDB.PDBParser(QUIET=True).get_structure(filename[0:4], filename)
        for chain in pdb.get_chains():
            io.set_structure(chain)
            io.save('./lanci/' + pdb.get_id() + chain.get_id() + '.pdb')

for files in os.listdir('./lanci/'):
    if files.endswith('.pdb'):
        os.system('gesamt * -o ' + files + '_out.pdb')            
        
