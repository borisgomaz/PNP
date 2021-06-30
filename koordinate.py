# coding: utf-8
import re
import pandas as pd
from Bio import PDB
f = open('4tta.tsv')
result = []
for n in f.readlines():
    r = re.match(r'(?P<frame>\d+)\s+(?P<bond_type>\w+)\s+(?P<chainID1>\w)\:(?P<AA1>\w+)\:(?P<AAID1>\d+)\:(?P<AAatom1>\w+)\s+(?P<chainID2>\w)\:(?P<AA2>\w+)\:(?P<AAID2>\d+)\:(?P<AAatom2>\w+)', n)
    rw = re.match(r'(?P<frame>\d+)\s+(?P<bond_type>\w+)\s+(?P<chainID1>\w)\:(?P<AA1>\w+)\:(?P<AAID1>\d+)\:(?P<AAatom1>\w+)\s+(?P<chainID2>\w)\:(?P<AA2>\w+)\:(?P<AAID2>\d+)\:(?P<AAatom2>\w+)\s+(?P<chainIDw>\w)\:(?P<water>\w+)\:(?P<waterID>\d+)\:(?P<water_atom>\w+)', n)
    rw2 = re.match(r'(?P<frame>\d+)\s+(?P<bond_type>\w+)\s+(?P<chainID1>\w)\:(?P<AA1>\w+)\:(?P<AAID1>\d+)\:(?P<AAatom1>\w+)\s+(?P<chainID2>\w)\:(?P<AA2>\w+)\:(?P<AAID2>\d+)\:(?P<AAatom2>\w+)\s+(?P<chainIDw>\w)\:(?P<water>\w+)\:(?P<waterID>\d+)\:(?P<water_atom>\w+)\s+(?P<chainIDw2>\w)\:(?P<water2>\w+)\:(?P<waterID2>\d+)\:(?P<water_atom2>\w+)', n)
    if rw2:
        result.append(rw2.groupdict())
    elif rw:
        result.append(rw.groupdict())
    elif r:    
        result.append(r.groupdict())
#df = pd.DataFrame(result)

parser = PDB.PDBParser()
pdb = '../4tta.pdb'
structure = parser.get_structure('4tta', pdb)
model = structure[0]
for res in result:
    chain1 = model[res.get('chainID1')]
    chain2 = model[res.get('chainID2')]
    residue1 = chain1[int(res.get('AAID1'))] 
    residue2 = chain2[int(res.get('AAID2'))]
    atom1 = residue1[res.get('AAatom1')] 
    atom2 = residue2[res.get('AAatom2')]
    
    distance = atom1 - atom2
    res.update({'distance': distance})
    
    if 'water' in res.keys():
        chain3 = model[res.get('chainIDw')]
        residue3 = chain3['W', int(res.get('waterID')), ' ']
        atom3 = residue3[res.get('water_atom')]
        
        distance_aw = atom1 - atom3
        distance_wa = atom3 - atom2
        res.update({'distance_atom1_water': distance_aw})
        res.update({'distance_water_atom2': distance_wa})
        
#Types of bonds: ["hp", "sb", "pc", "ps", "ts", "vdw", "hb", "wb", "wb2"]
#For hb bonds: if m['bond_type'][:2] == 'hb':
short = []
for m in result:
    if m['bond_type'] == 'wb':
        short.append(m)
df = pd.DataFrame(short)
#For maximal rows dislay:
#pd.set_option('display.max_rows', df.shape[0]+1)
print(df)        
    
