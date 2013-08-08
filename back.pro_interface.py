#! /usr/bin/env python
#interface for pro_angle.py
#by Marcio von Muhlen, summer 2004

# USAGE:
#
#
# from directory containing *.ent files
# ]$ python pro_interface [res seq] [seq index] [output name] //cont
# [chi on/off] [input file] [output_format
#
# [res seq] is the residue sequence to look for, format 'GLY PRO ARG ..'
# if desire to look at all residues, use 'all'
# [seq index] specifies which residue to return angles for, e.g. '1' for PRO
# above (python indexes by 0)
# [output name] is name of file to write output to
# [chi on/off] specifies whether to calculate only phi/psi angles
# or also chi1 angles in output. Use n for phi/psi, y for chi
#
# DEFAULTS: chi defaults to off, output name to angle_out.txt, seq index to 0,
# res seq to 'PRO'


def compile_angles(res = 'all', res_angle = 0, chain = 0,
                   file_out = 'angle_out.txt', chi_ = 'n',
                   input_file = 'dir', output_format = 'std'):
    """
    res: residue sequence
    
    res_angle: index of the residue in that sequence to look at

    chain: chain in PDB file
    
    file_out: filename of the output file
    
    chi:  look for phi psi only ('n') or also chi angles ('y')
    
    input_file: default is to look at all .ent and pdb iles in this directory

    output_format: default is gnuplot data format, 'alt' is labeled residues
    (only makes sense for res = 'all', else ignored)
    """
    import commands
    from pro_angle import pro_angle
    if input_file == 'dir':
        file_list = commands.getoutput('ls *.ent').split()
    else:
        file_list = [input_file]
    file_out = open(file_out, 'w')
    if chi_ == 'n' and res != 'all':
        header = '# phi, psi angles for residue #'  + str(res_angle)
    elif res == 'all' and chi_ == 'n':
        header = '#phi, psi angles for all residues'
    elif res == 'all' and chi != 'n':
        header = '#phi, psi, chia and chib  angles for all residues'        
    else:
        header = '# phi, psi, chia, chib angles for residue #' + str(res_angle)
    if res_angle > (len(res.split())+1):
        raise IndexError   
    file_out.write(header + ' in sequence: ' +
                   res + '\n#from file(s): ' + input_file +'\n')
    count = 0
    if res == 'all' and output_format != 'std':
        detailed_output(file_out, chi_, file_list, output_format)
    else:   
        try:
            if res_angle > (len(res.split())-1):
                raise IndexError        
            for pdb in file_list:
                print "processing: " + str(pdb)
                angle_list=pro_angle(pdb, res, chain, 0, res_angle, chi_)
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
            print "\nINPUT ERROR: check to make sure .ent or .pdb files are present, res_angle < len(res_seq)\n"

def detailed_output(file_out, chi_, file_list, output_format):
    aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
           'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
           'ARG', 'HIS']
    count = 0
    for pdb in file_list:
        print "processing: " + str(pdb)
        for aa in aas:
            angle_list=pro_angle(pdb, aa, 0, 0, 0, chi_)
            if not angle_list: #no angles found for this residue               
                file_out.write('#Residue: ' + aa + ' Angles: NULL\n')
            else:
                #print res id, print values
                file_out.write('#Residue: ' + aa + ' Angles: \n')
                for tuple in angle_list:
                    if len(tuple) > 2:
                        file_out.write(str(tuple[0]).zfill(6) + '  '
                                       +str(tuple[1]).zfill(6)+ '  '
                                       +str(tuple[2]).zfill(6)+ '  '
                                       +str(tuple[3]).zfill(6)+' ')
                    else:
                        file_out.write(str(tuple[0]).zfill(6)+ '  '
                                       +str(tuple[1]).zfill(6)+'\n')
                count += 1
    file_out.write('# of  residues: ' + str(count)+'\n')
    file_out.close()
    

##main

import sys
import string

if len(sys.argv) > 1:
    res = sys.argv[1]
else:
    res = 'PRO'
if len(sys.argv) > 2:
    ang = string.atoi(sys.argv[2])
else:
    ang = 0
if len(sys.argv) > 3:
    chain = sys.argv[3]
else:
    chain = 'default' 
if len(sys.argv) > 4:
    file_out = sys.argv[4]
else:
    file_out = "angle_out.txt"
if len(sys.argv) > 5:
    chi = sys.argv[5]
else:
    chi = 'n'
if len(sys.argv) > 6:
    input = sys.argv[6]
else:
    input = 'dir'
if len(sys.argv) > 7:
    output = sys.argv[7]
else:
    output = 'std'

print 'parsed \nres: ' + res + '\nangle_index: ' + str(ang) \
      + '\nchain: ' + chain + '\nfile_out: ' + file_out + '\nchi? ' \
      + chi + '\ninput_file? ' + input

compile_angles(res, ang, chain,  file_out, chi, input, output)



    
