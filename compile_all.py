#! /usr/bin/env python
#script to get dihedral angle data from pdb files into one file

def detailed_output(file_out, chi_, file_list, chain = 0):
    #added last parameter, chain, on 7/28
    """to be used by compile angles only
    from pdbs in a file list, puts phi psi (chi?) angles in file_out
    """

    os.chdir('/home/marciovm/proteins/')
    from pro_angle import pro_angle
    os.chdir('/home/marciovm/proteins/dist')
    
    aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
           'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
           'ARG', 'HIS']
    count = 0
    for pdb in file_list:
        print "processing: " + str(pdb)
        for aa in aas:
            angle_list=pro_angle(pdb, aa, chain, 0, 0, chi_)
            if not angle_list: #no angles found for this residue               
                file_out.write('\n#Residue ' + aa + ' angles: NULL\n')
            else:
                file_out.write('\n#Residue ' + aa + ' angles: \n')
                for tuple in angle_list:
                    if len(tuple) > 2:
                        file_out.write(str(tuple[0]).zfill(6) + '  '
                                       +str(tuple[1]).zfill(6)+ '  '
                                       +str(tuple[2]).zfill(6)+ '  '
                                       +str(tuple[3]).zfill(6)+'\n')
                    else:
                        file_out.write(str(tuple[0]).zfill(6)+ '  '
                                       +str(tuple[1]).zfill(6)+'\n')
                    count += 1
    file_out.write('\n# of  residues: ' + str(count)+'\n')
    file_out.close()



import os
import commands
import sys
import string

os.chdir('/home/marciovm/proteins/dist')

list = commands.getoutput('ls *.ent').split()

start = sys.argv[1]
end = sys.argv[2]

count = 0
for a in range(string.atoi(start), string.atoi(end)):
    print count
    count += 1
    pdb = list[a]
    filename = '/home/marciovm/proteins/dist/data/' + pdb + '.dat'
    file_out = open(filename, 'w')
    if pdb[7] == '_':
       detailed_output(file_out, 'y', [pdb], ' ')
    else:
       detailed_output(file_out, 'y', [pdb], pdb[7])

   


    
