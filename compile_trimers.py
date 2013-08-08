#! /usr/bin/env python
#compile trimer phi psi data
#(that means 8k rama plots)

#chain number fetch is weird!!

#! /usr/bin/env python
#script to get dihedral angle data from pdb files into trimer data files
# main

def compile_angles(res = 'all', res_angle = 0, chain = 0,
                   file_out = 'angle_out.txt', chi_ = 'n',
                   input_file = 'dir', output_format = 'std',
                   sec_str = 'all', reverse = 1):
    """
    """
    import sys
    import commands
    os.chdir('/home/marciovm/proteins')
    from pro_angle import pro_angle
    os.chdir('/home/marciovm/proteins/dist')
    file_list = [input_file]
    file_out = open(file_out, 'a')
    count = 0
    if res == 'all' and output_format != 'std':
        pass
    else:   
        try:
            if res_angle > (len(res.split())-1):
                raise IndexError        
            for pdb in file_list:
                angle_list=pro_angle(pdb, res, chain, 0, res_angle,
                                     chi_, sec_str, reverse)
                for tuple in angle_list:
                    if tuple:
                        if len(tuple) > 2:
                            file_out.write(str(tuple[0]).zfill(6) + '  '
                                           +str(tuple[1]).zfill(6)+ '  '
                                           +str(tuple[2]).zfill(6)+ '  '
                                           +str(tuple[3]).zfill(6)+'\n')
                        else:
                            file_out.write(str(tuple[0]).zfill(6)+ '  '
                                           +str(tuple[1]).zfill(6)+'\n')
                        count += 1
                    else:
                        pass
            file_out.write('#occurences of sequence: ' + str(count)+'\n')
            file_out.close()
        except IndexError:
            errormsg =  '\nINPUT ERROR: check to make sure .ent or .pdb files are present, res_angle < len(res_seq)\n'
            print errormsg
            

#main
import os
import commands
import sys
import string

os.chdir('/home/marciovm/proteins/dist')
list = commands.getoutput('ls *.ent').split() #list of pdb files w/ chain names
start = sys.argv[1]
end = sys.argv[2]
count = 0

aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
       'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
       'ARG', 'HIS']

for pre_index in range(string.atoi(start), string.atoi(end)):
    pre = aas[pre_index]
    for cur in aas:
        for post in aas:
            sequence = pre + ' ' + cur + ' ' + post
            filename = '/home/marciovm/proteins/trimers/' + string.lower(cur) + '/' + pre + "_" + cur + "_" + post + '.dat'
            print '\nlooking at ' + sequence
            count += 1
            for pdb in list:
                file = open(filename, 'w')
                header = '#phi, psi, chia, chib angles for m-res in ' + sequence+ '\n'
                file.write(header)
                file.close()
                
                if pdb[7] == '_':
                    compile_angles(sequence, 1, ' ', filename, 'y', pdb)
                else:
                    compile_angles(sequence, 1, pdb[7], filename, 'y', pdb)

