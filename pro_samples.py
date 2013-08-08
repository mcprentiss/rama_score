#test file for pro_angle

aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
       'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
       'ARG', 'HIS']

p = aas[14]

from pro_interface import *

for res in aas:
    filename = 'PRO-' + res + '.txt'
    seq = 'PRO ' + res
    print "\namino acid sequence: " + seq + '\n'
    compile_angles(seq, 0, filename, 'y')


#interesting: pro-pro, pro-asn

