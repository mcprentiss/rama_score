#! /usr/bin/env python
#A collection of functions that calculate the m-score for a given PDB file
#uses score_table files in working directory
#by Marcio von Muhlen, summer 2004

def score(PDBfile):
    """
    Calculates the m-score for a given PDB file

    arguments:
    
    PDBfile - the PDB file to score

    hidden arguments:

    aas.scr, pro.scr, gly.scr - the scoring tables
    need to be present in working directory
    
    """
    from pro_angle import find_residue
    from Bio.PDB.PDBParser import PDBParser
    from pro_length import length
    
    (aas, gly, pro) = load_scores() ##define global tables
    score = 0 #initialize    
    pars = PDBParser(PERMISSIVE = 1)
    struct = pars.get_structure(PDBfile.rstrip('.pdb'), PDBfile)
    model = struct.child_list[0]
    chain = model.child_list[0]
    pro_list = find_residue(chain, 'PRO')
    gly_list = find_residue(chain, 'GLY')
    aas_list = range(chain.child_list[1].id[1],
                     chain.child_list[len(chain)-1].id[1])
    #need to remove pro/gly indices in first/last position
    if pro_list.count(1) > 0:        
        pro_list.remove(1)
    if pro_list.count(len(chain)-1) > 0:
        pro_list.remove(len(chain)-1)
    if gly_list.count(1) > 0:
        gly_list.remove(1)
    if gly_list.count(len(chain)-1) > 0:
        gly_list.remove(len(chain)-1)   
    try:
        for index in pro_list:       
            aas_list.remove(index) #remove pros from aas_list
        for index in gly_list:
            aas_list.remove(index) #remove glys from aas_list
    except ValueError:
        print 'incosistency in PDB file - will return score = 0' 
        return 0
    else:
        proscore = score_help(chain, pro_list, pro)
        glyscore = score_help(chain, gly_list, gly)
        aasscore = score_help(chain, aas_list, aas)
        score = proscore+glyscore+aasscore
        size=length(chain)
        try:
            score = (score/size)*1000 #normalize score
            return score
        except ZeroDivisionError:
            print "calculated protein length 0 -> returning score 0"
            score = 0
            return score
  

def score_help(chain, res_list, score_table):
    """
    Helper function

    Compares angles of given residue list to score_table

    Returns a score

    Arguments:

    chain - name of the chain in the PDB model to examine

    res_list - list of residue indices

    score_table - table to look up scores in (meant to be residue specific)
    """
    from pro_angle import calc_dihedral
    from math import floor
    
    score = float(0)
    for res in res_list:
        try:
            (phi,psi) = calc_dihedral(chain, res) 
            indx = int(floor(phi/10)+18)
            indy = int(floor(psi/10)+18)
            temp = float(score_table[indy][indx])
            score = score + temp
        except ValueError:
            pass
#            print "ValueError: asked for score of non-scorable residue"
    return score

def load_scores():
    """
    Helper function, loads score tables from files

    Returns score tables    
    """
    from copy import copy
    from string import atof
    aas = open('aas.scr')
    pro = open('pro.scr')
    gly = open('gly.scr')

    aasline = aas.readline().split()
    proline = pro.readline().split()
    glyline = gly.readline().split()
    
    probx = [0 for i in xrange(36)] #this will be x index
    proby = [0 for i in xrange(36)] #this will be y index    

    for row_counter in range(36):
        for column_counter in range(36):
            probx[column_counter] = atof(aasline[column_counter])
        aasline = aas.readline().split()
        proby[row_counter] = copy(probx)
    aas = copy(proby)

    probx = [0 for i in xrange(36)]
    proby = [0 for i in xrange(36)]
    for row_counter in range(36):
        for column_counter in range(36):
            probx[column_counter] = atof(proline[column_counter])
        proline = pro.readline().split()
        proby[row_counter] = copy(probx)
    pro = copy(proby)

    probx = [0 for i in xrange(36)]
    proby = [0 for i in xrange(36)]
    for row_counter in range(36):
        for column_counter in range(36):
            probx[column_counter] = atof(glyline[column_counter])
        glyline = gly.readline().split()
        proby[row_counter] = copy(probx)
    gly = copy(proby) 
    return (aas, gly, pro)
        
#main

#import sys
#pdb = sys.argv[1]
#a=score(pdb)
#print a
    

