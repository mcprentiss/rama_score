from Bio.PDB.PDBParser import PDBParser
p = PDBParser(PERMISSIVE = 1)
struct = p.get_structure('t', 'egp')
model = struct.child_list[0]
chain = model.child_list[0]
