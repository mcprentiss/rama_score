#! /usr/bin/env python
#A command line interface for pro_angle.py
#by Marcio von Muhlen, summer 2004


def compile_angles(res = 'all', res_angle = 0, chain = 0,
                   file_out = 'angle_out.txt', chi_ = 'n',
                   input_file = 'dir', output_format = 'std',
                   sec_str = 'all', reverse = 1, start_index = 0, end_index =
                   10000):
    """

    A command line interface function for pro_angle

    Outputs calculated angles to a specified output file
    
    Arguments:
    
    res - residue sequence
    'all' means every residue
    
    res_angle - index of the residue in sequence to examine

    chain - name of chain in PDB file to examine
    
    file_out - filename of the output file
    
    chi - on/off switch for chi angle output (use y/n)
        
    input_file - filename of the input PDB file
    If unspecified, look at all *.ent files in working directory

    output_format - syntax of output
    'alt' outputs angles categorized by residue 
    'std' outputs unlabeled angles

    'alt' only makes sense if res = 'all'

    sec_str - limits search to DSSP secondary structure identifier

    reverse - limits search to all non-sec_str residues
    
    DSSP codes:
      - H        Alpha helix (4-12)
      - B        Isolated beta-bridge residue
      - E        Strand
      - G        3-10 helix
      - I        pi helix
      - T        Turn
      - S        Bend
      - -        None

    start_index - limits search to residues w/ chain numbers above this one

    end_index - limits search to residues w/ chain numbers below this one
    """
    import commands
    from pro_angle import pro_angle

    sec = '\n#sec_str: ' + str(sec_str) + ' excluded? ' + str(reverse)
    
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
                   res + '\n#from file(s): ' + input_file + sec +'\n')
    count = 0
    if res == 'all' and output_format != 'std':
        detailed_output(file_out, chi_, file_list)
    else:   
        try:
            if res_angle > (len(res.split())-1):
                raise IndexError        
            for pdb in file_list:
                print "processing: " + str(pdb)
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
            print "\nINPUT ERROR: check to make sure .ent or .pdb files are present, res_angle < len(res_seq)\n"

def detailed_output(file_out, chi_, file_list):
    """
    Helper Function    
    best left to usage by compile_angles function

    Arguments:
    file_out - file_out - filename of the output file

    file_list - list of PDB filenames

    chi_ - on/off switch for chi angle output (use y/n)
    """

    from pro_angle import pro_angle
    aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
           'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
           'ARG', 'HIS']
    count = 0
    for pdb in file_list:
        print "processing: " + str(pdb)
        for aa in aas:
            angle_list=pro_angle(pdb, aa, 0, 0, 0, chi_)
            if not angle_list: #no angles found for this residue               
                file_out.write('\n#Residue ' + aa + ' angles: NULL\n')
            else:
                file_out.write('\n#Residue ' + aa + ' angles: \n')
                for tuple in angle_list:
                    if len(tuple) > 2:
                        file_out.write('(' + str(tuple[0]).zfill(6) + '  '
                                       +str(tuple[1]).zfill(6)+ '  '
                                       +str(tuple[2]).zfill(6)+ '  '
                                       +str(tuple[3]).zfill(6)+') ')
                    else:
                        file_out.write('(' + str(tuple[0]).zfill(6)+ '  '
                                       +str(tuple[1]).zfill(6)+') ')
                    count += 1
    file_out.write('\n# of  residues: ' + str(count)+'\n')
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
if len(sys.argv) > 8:
    sec_str = sys.argv[8]
else:
    sec_str = 'all'
if len(sys.argv) > 9:
    rev = sys.argv[9]
else:
    rev = 1
    
if sec_str != 'all':
    print 'parsed \nres: ' + res + '\nangle_index: ' + str(ang) \
      + '\nchain: ' + chain + '\nfile_out: ' + file_out + '\nchi? ' \
      + chi + '\ninput_file? ' + input + '\nsec_str? ' + sec_str + \
      '\nreverse? ' + str(rev)
else:
    print 'parsed \nres: ' + res + '\nangle_index: ' + str(ang) \
      + '\nchain: ' + chain + '\nfile_out: ' + file_out + '\nchi? ' \
      + chi + '\ninput_file? ' + input + '\nsec_str? all'
    
compile_angles(res, ang, chain, file_out, chi, input, output, sec_str, rev)


    
