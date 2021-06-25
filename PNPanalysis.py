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

#To ignore warnings: QUIET=True
parser = PDB.PDBParser(QUIET=True)
pdb = '../4tta.pdb'
structure = parser.get_structure('4tta', pdb)
for res in result:
    atom1 = structure[0][res.get('chainID1')][int(res.get('AAID1'))][res.get('AAatom1')] 
    atom2 = structure[0][res.get('chainID2')][int(res.get('AAID2'))][res.get('AAatom2')]
        
    if 'water2' in res.keys():
        atom_w1 = structure[0][res.get('chainIDw')]['W', int(res.get('waterID')), ' '][res.get('water_atom')]
        atom_w2 = structure[0][res.get('chainIDw2')]['W', int(res.get('waterID2')), ' '][res.get('water_atom2')]
          
        distance_a1w1 = atom1 - atom_w1
        distance_w1w2 = atom_w1 - atom_w2
        distance_w2a2 = atom_w2 - atom2
        res.update({'distance_atom1_water1': distance_a1w1})
        res.update({'distance_water1_water2': distance_w1w2})
        res.update({'distance_water2_atom2': distance_w2a2})
    
    elif 'water' in res.keys():
        atom_w = structure[0][res.get('chainIDw')]['W', int(res.get('waterID')), ' '][res.get('water_atom')]
      
        distance_aw = atom1 - atom_w
        distance_wa = atom_w - atom2
        res.update({'distance_atom1_water': distance_aw})
        res.update({'distance_water_atom2': distance_wa})
    
    else:
        distance = atom1 - atom2
        res.update({'distance': distance})

#Types of bonds: ["hp", "sb", "pc", "ps", "ts", "vdw", "hb", "wb", "wb2"]
#For hb bonds: if m['bond_type'][:2] == 'hb':
short = []
for m in result:
    if m['bond_type'] == 'sb':
        short.append(m)
df = pd.DataFrame(short)
#For maximal rows dislay:
#pd.set_option('display.max_rows', df.shape[0]+1)
print(df)        
        
