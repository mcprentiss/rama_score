#! /usr/bin/env python
#compile trimer phi psi data
#(that means 8k rama plots)

def compile_trimers():

    import os
    import commands
    import sys
    import string
    from Bio.PDB.PDBParser import PDBParser
    from pro_length import length    
    
    os.chdir('/home/marciovm/proteins/bigdist') ##choose directory
    
    list = commands.getoutput('ls *.ent').split() #list of pdb files w/ chain names
    
    start = sys.argv[1]
    end = sys.argv[2]
    
    aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
           'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
           'ARG', 'HIS']
    
    for pdb_index in range(string.atoi(start), string.atoi(end)):
        pdb = list[pdb_index]

    
        pars = PDBParser(PERMISSIVE = 1)
        struct = pars.get_structure(pdb, pdb)
        model = struct.child_list[0]
        if pdb[4] == '.':
            chain_name = ' '
        else:
            chain_name = string.upper(pdb[4])            
        chain = model.child_dict[chain_name]

        print 'processing ' + pdb + ' index: ' + str(pdb_index)

        ##find first real residue
        index = 0
        res = 'first'
        while aas.count(res) == 0 and index < 10000:
            try:
                res = chain.child_list[index].resname
            except KeyError:
                pass #nothing here
            index += 1
        res = chain.child_list[index].resname ##should be 2nd real residue
        ##find last real residue
        final_res_index = length(chain) + index - 3 ##2nd to last real residue
        ##look at everything in between
        for aa_index in range(index, final_res_index):
            aacur = chain.child_list[aa_index].resname
            aapre = chain.child_list[aa_index - 1].resname
            aapost= chain.child_list[aa_index + 1].resname
            seq = aapre + '_' + aacur + "_" + aapost + '.dat'
            filename = ('/home/marciovm/proteins/bdtrimers/' + string.lower(aacur) + '/' + seq)
            try:
                file_out = open(filename, 'a')
            except IOError:
                print "IOError line 63, nostandard aa name in PDB"
            else:
                tuple = calc_all_dihedrals(chain, aa_index)
                if tuple:
                    if len(tuple) > 2:
                        file_out.write(str(tuple[0]).zfill(6) + '  '
                                       +str(tuple[1]).zfill(6)+ '  '
                                       +str(tuple[2]).zfill(6)+ '  '
                                       +str(tuple[3]).zfill(6)+'\n')
                    else:
                        file_out.write(str(tuple[0]).zfill(6)+ '  '
                                       +str(tuple[1]).zfill(6)+'\n')
                file_out.close()
                

                
def calc_all_dihedrals(chain, child_res_id):
    """
    calculates the dihedral angles (phi, psi, chia, chib) for residue of
    index child_res_id in chain
    
    returns a tuple of form (phi, psi, chia, chib), if it exists
    
    where ambiguity in chi angle definition exists:
    chia is in reference to the longer side chain or the heavier atom
    chib to the shorter
    
    if no ambiguity, chia=chib
    
    if residue is a GLY or ALA return only (phi, psi)
    """

    from math import pi
    from Bio import PDB
    try:
        residue = chain.child_list[child_res_id]
        name = residue.get_resname()
        if name == 'GLY' or name == 'ALA':
            dih = calc_dihedral(chain, child_res_id)
            phi = dih[0]
            psi = dih[1]
            return (phi, psi, 0, 0)          
        
    except KeyError:
        print "key error line 103"
        return () # no dihedral angles for corner residues or non-a.a.'residues'
    except IndexError:
        print 'IndexError line 106, probable cause: irregular PDB file'
        return ()     
    try:
        CP = chain.child_list[(child_res_id-1)]['C'].get_vector()
        N = chain.child_list[child_res_id]['N'].get_vector()
        CA = chain.child_list[child_res_id]['CA'].get_vector()
        C = chain.child_list[child_res_id]['C'].get_vector()
        NA = chain.child_list[(child_res_id+1)]['N'].get_vector()
        CB = chain.child_list[child_res_id]['CB'].get_vector()
        fourth_chi_atom = chain.child_list[child_res_id].child_list[5].get_vector()
        if name == 'VAL' or name == 'ILE' or name == 'THR':
            alt_fourth_chi_atom = chain.child_list[child_res_id].child_list[6].get_vector()
        else:
            alt_fourth_chi_atom = fourth_chi_atom
    except KeyError:
        print 'KeyError line 119'
        return () # no dihedral angles for corner residues or non-a.a.'residues'
    except IndexError:
        print 'IndexError line 122, probable cause: irregular PDB file'
        return ()
    else:
        try:
            phi = PDB.calc_dihedral(CP, N, CA, C)*-180/pi
            psi = PDB.calc_dihedral(N, CA, C, NA)*-180/pi
            chia= PDB.calc_dihedral(C, CA, CB, fourth_chi_atom)*-180/pi
            chib= PDB.calc_dihedral(C, CA, CB, alt_fourth_chi_atom)*-180/pi
            return (phi, psi, chia, chib)
        except ZeroDivisionError:
            return ()



def calc_dihedral(chain, child_res_id):
    """
    calculates the dihedral angles (phi, psi) for residue of index
    child_res_id in chain
    
    returns a tuple of form (phi, psi), if it exists
    
    """
    from math import pi
    from Bio import PDB
    
    try:
        CP = chain.child_list[(child_res_id-1)]['C'].get_vector()
        N = chain.child_list[child_res_id]['N'].get_vector()
        CA = chain.child_list[child_res_id]['CA'].get_vector()
        C = chain.child_list[child_res_id]['C'].get_vector()
        NA = chain.child_list[(child_res_id+1)]['N'].get_vector()
    except KeyError:
        return () # no dihedral angles for corner residues or non-a.a.'residues'
    else:
        try:
            phi = PDB.calc_dihedral(CP, N, CA, C)*-180/pi
            psi = PDB.calc_dihedral(N, CA, C, NA)*-180/pi
            return (phi, psi)
        except ZeroDivisionError:
            return ()


#main
compile_trimers()
 
