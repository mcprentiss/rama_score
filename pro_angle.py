
#A collection of functions which utilize biopython to analyze torsion angles
#in PDB files
#by Marcio von Muhlen, summer 2004


def pro_angle(PDBfile='egp', res_name='GLY', chain_name = 'default',
              model_num = 0, res_out_angle = 0, chi = 'y',
              sec_str = 'all', reverse = 0, start_index = 0, end_index = 9999):
    """
    
    Calculates torsion angles in input pdb file
     
    Arguments:

    PDBfile - the input file to examine

    res_name - the sequence of residues to look for, i.e. 'PRO ARG'
    
    chain_name - name of the chain in the PDB model to examine

    model_num - name of the model in the PDB file to examine

    res_out_angle - index of the residue in the sequence to examine

    chi - on/off switch for chi angle output (use y/n)

    sec_str - limits search to DPSS secondary structure identifier

    reverse - limits search to all non-sec_str residues

    start_index - limits search to residues w/ chain numbers above this one

    end_index - limits search to residues w/ chain numbers below this one
    """
    
    from Bio.PDB.PDBParser import PDBParser
    
    pars = PDBParser(PERMISSIVE = 1)
    struct = pars.get_structure(PDBfile.rstrip('.pdb'), PDBfile)
    model = struct.child_list[model_num]
    if chain_name == 'default' or chain_name == 0 or chain_name == '0':
        chain = model.child_list[0]
    else:    
        chain = model.child_dict.get(chain_name)
       
    res_list = res_name.split()
    if res_list[0] != 'all':       
        seq_index = find_residue(chain, res_list[0]) #narrow to indices of first res
        res_list.remove(res_list[0])
        seq_index = find_sequence(chain, res_list, seq_index) #find desired seqs    
    else: #'all' case


        if end_index > len(chain):
            end_index = len(chain)
        
        ####### changed
        seq_index = range(start_index, end_index)
        ####### changed

    if sec_str != 'all': #run filter to eliminate residues with incorrect sec_str
        if reverse == 0 or reverse == '0':
            seq_index = sec_str_filter(PDBfile, model, chain, seq_index, sec_str)
        else:
            seq_index = sec_str_rev_filter(PDBfile, model, chain, seq_index, sec_str)
    dihedrals = []
    if chi == 'n':
        try:
            for res in seq_index:
                angles = calc_dihedral(chain, res + res_out_angle)
                if angles:
                    phi = round(angles[0])
                    psi = round(angles[1])
                    dihedrals.append((phi, psi))
                else:
                    #dihedrals.append(()) #don't want blanks
                    pass
        except IndexError:
            pass
    else:
        for res in seq_index:              
            try:
                angles = calc_all_dihedrals(chain, res + res_out_angle)
                if angles:
                    phi = round(angles[0])                    
                    psi = round(angles[1])                    
                    chia = angles[2]
                    if chia == 0:
                        dihedrals.append((phi, psi))                                            
                    else:
                        chia = round(angles[2])
                        chib = round(angles[3])
                        dihedrals.append((phi, psi, chia, chib))
                else:                    
                    #dihedrals.append(())
                    pass
            except IndexError:
                #print "IndexError 86 - Tried to calculate dihedral of terminal residue"
                pass
           # pass #first and last residues will cause errors, expected
    
    return dihedrals
   

def find_sequence(PDBchain, res_list, seq_index, count = 0):
    """       
    Returns list of start indexes for residue sequence res_list, recursive

    looks at chain.child_list indices, NOT chain indices!
    
    if no instance found returns 0
    
    seq_index -  list of locations of first residue
    
    res_list  -  list of sequence starting w/ 2nd residue 

    """
    if not res_list: #base case, nothing left on res_list
        return seq_index
    count += 1
    r = 0
    while r < len(seq_index):
        try:
            if res_list[0] == 'XXX':
                r += 1
            elif PDBchain.child_list[seq_index[r]+count].get_resname() != res_list[0]:
                seq_index.remove(seq_index[r]) ##remove this search path   
            else:
                r += 1
        except KeyError: #searching for a non-existant residue
            seq_index.remove(seq_index[r])          
    res_list.remove(res_list[0]) #this res found already, can go away          
    return find_sequence(PDBchain, res_list, seq_index, count = count)
            
                
def find_residue(chain, res_name):
    """
    Looks for occurences of residue res_name in chain
    returns a list with id numbers of matches, if any
    
    looks at chain.child_list indices, NOT chain indices!
    """    
    match_list = []
    iter = chain.get_iterator()
    current = iter.next()
    count = 0
    while current:
        if current.get_resname() == res_name:
            match_list.append(count)
        try:
            current = iter.next()
            count += 1
        except StopIteration:
            current = 0
    return match_list


def calc_dihedral(chain, res_id):
    """
    calculates the dihedral angles (phi, psi) for residue of index
    res_id in chain
    
    returns a tuple of form (phi, psi), if it exists    
    """
    
    from math import pi
    from Bio import PDB
    
    try:
        CP = chain.child_list[(res_id-1)]['C'].get_vector()
        N = chain.child_list[res_id]['N'].get_vector()
        CA = chain.child_list[res_id]['CA'].get_vector()
        C = chain.child_list[res_id]['C'].get_vector()
        NA = chain.child_list[(res_id+1)]['N'].get_vector()       
    except KeyError:
        return () # no dihedral angles for corner residues or non-a.a.'residues'
    except IndexError:
        return () # no dihedral angles for corner residues or non-a.a.'residues'
    else:
        phi = PDB.calc_dihedral(CP, N, CA, C)*-180/pi
        psi = PDB.calc_dihedral(N, CA, C, NA)*-180/pi
        return (phi, psi)


def calc_all_dihedrals(chain, res_id):
    """
    calculates the dihedral angles (phi, psi, chia, chib) for residue of
    index res_id in chain
    
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
        residue = chain.child_list[res_id]
        name = residue.get_resname()
        if name == 'GLY' or name == 'ALA':
            dih = calc_dihedral(chain, res_id)
            phi = dih[0]
            psi = dih[1]
            return (phi, psi, 0, 0)          
        
    except KeyError:
        print "key error line 196"
        return () # no dihedral angles for corner residues or non-a.a.'residues'    
    try:
        CP = chain.child_list[(res_id-1)]['C'].get_vector()
        N = chain.child_list[res_id]['N'].get_vector()
        CA = chain.child_list[res_id]['CA'].get_vector()
        C = chain.child_list[res_id]['C'].get_vector()
        NA = chain.child_list[(res_id+1)]['N'].get_vector()
        CB = chain.child_list[res_id]['CB'].get_vector()
        fourth_chi_atom = chain.child_list[res_id].child_list[5].get_vector()
        
        if name == 'VAL' or name == 'ILE' or name == 'THR':
            alt_fourth_chi_atom = chain.child_list[res_id].child_list[6].get_vector()
        else:
            alt_fourth_chi_atom = fourth_chi_atom
    except KeyError:
        print 'key error line 211'
        return () # no dihedral angles for corner residues or non-a.a.'residues'
    else:
        phi = PDB.calc_dihedral(CP, N, CA, C)*-180/pi
        psi = PDB.calc_dihedral(N, CA, C, NA)*-180/pi
        chia= PDB.calc_dihedral(C, CA, CB, fourth_chi_atom)*-180/pi
        chib= PDB.calc_dihedral(C, CA, CB, alt_fourth_chi_atom)*-180/pi
    return (phi, psi, chia, chib)


def sec_str_filter(PDBfile, model, chain, seq_index, sec_str):
    """
    Helper function

    A filter for secondary structures
    
    Checks residues in seq_index of chain, removes those w/non-specd sec_str

    Returns filtered seq_index
    """
    from Bio.PDB import DSSP
    dssp = DSSP(model, PDBfile) ##need to hide this output
    num = 0
    while num < len(seq_index):
        try:
            res = chain.child_list[seq_index[num]]
            sec = dssp.__getitem__(res)[0]    
            if sec != sec_str:
                seq_index.remove(seq_index[num])
            else:            
                num += 1
        except KeyError: #not an amino acid
            num += 1
    return seq_index

def sec_str_rev_filter(PDBfile, model, chain, seq_index, sec_str):
    """
    Helper function

    A filter for secondary structures.

    Reverse filter of sec_str_filter
    
    Checks residues in seq_index of chain, removes those w/ specified sec_str

    Returns filtered seq_index    
    """
    from Bio.PDB import DSSP
    dssp = DSSP(model, PDBfile) ##need to hide this output
    num = 0
    while num < len(seq_index):
        try:
            res = chain.child_list[seq_index[num]]
            sec = dssp.__getitem__(res)[0]    
            if sec == sec_str:
                seq_index.remove(seq_index[num])
            else:            
                num += 1
        except KeyError: #not an amino acid
            num += 1
    return seq_index
    
