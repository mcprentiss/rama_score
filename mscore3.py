#! /usr/bin/env python
#A collection of functions that calculate the m-score for a given PDB file
#uses score_table files in working directory
#by Marcio von Muhlen, summer 2004

#update 2: uses trimer score tables (N=8k) instead of pro, gly, other
#update 3: possible to use mutation list 


def score(PDBfile, start_index = 1, end_index = 9999, smoothed = 'y', mutated_list = ' '):
    """
    Calculates the m-score for a given PDB file

    arguments:
    
    PDBfile - the PDB file to score

    start_index - child_list stat index of residues to look at

    end_index - child_list end index
    
    hidden arguments:

    aas.scr, pro.scr, gly.scr - the scoring tables
    need to be present in working directory
    
    """
    from pro_angle import find_residue
    from Bio.PDB.PDBParser import PDBParser
    from pro_length import length
    import os
    import string

    score = 0 #initialize    
    pars = PDBParser(PERMISSIVE = 1)
    struct = pars.get_structure(PDBfile.rstrip('.pdb'), PDBfile)
    model = struct.child_list[0]
    chain = model.child_list[0]

    score = float(0)
    size=length(chain)

    if end_index > (size - 2) :
        end_index = size - 2

    #non-mutated score
    if mutated_list == ' ':        
        for res_index in range(start_index, end_index): #not first or last res
            res = chain.child_list[res_index]
            cur = res.resname
            pre = chain.child_list[res_index-1].resname
            pos = chain.child_list[res_index+1].resname
            if smoothed == 'y':
                filename = pre + '_' + cur + '_' + pos + '.smooth.scr'
            else:
                filename = pre + '_' + cur + '_' + pos + '.scr'        
            table_file = '/home/marciovm/proteins/bdtrimers/' + \
            string.lower(cur) + '/' + filename
            chain_index = chain.child_list[res_index].id[1]            

#           print 'loading table file: ' + table_file ### debugging
            
            table = load_scores(table_file)
            if table != 0:
                new = score_help(chain, chain_index, table)
            else:
                new = 0
            score = score + new
    else:
        mutated_form = mutated_list.split()       
        for res_index in range(start_index, end_index):
            res = chain.child_list[res_index]
            mutated_index = res_index - start_index            
            cur = mutated_form[res_index-start_index]
            if mutated_index == 0:
                pre = chain.child_list[res_index-1].resname
            else:
                pre = mutated_form[mutated_index-1]
            if mutated_index == end_index - start_index - 1:
                pos = chain.child_list[res_index+1].resname
            else:
                pos = mutated_form[mutated_index+1]
            if smoothed == 'y':
                filename = pre + '_' + cur + '_' + pos + '.smooth.scr'
            else:
                filename = pre + '_' + cur + '_' + pos + '.scr'
                
            table_file = '/home/marciovm/proteins/bdtrimers/' + \
            string.lower(cur) + '/' + filename
            chain_index = chain.child_list[res_index].id[1]

            # print 'loading table file: ' + table_file ##debugging
            
            table = load_scores(table_file)
            if table != 0:
                new = score_help(chain, chain_index, table)
            else:
                new = 0
            score = score + new            
    try:
        score = (score/size)*1000 #normalize score
        return score
    except ZeroDivisionError:
        print "calculated protein length 0 -> returning score 0"
        score = 0
        return score
  

def score_help(chain, res_id, score_table):
    """
    Helper function

    Computes mscore for middle residue in given trimer

    Returns a score

    Arguments:

    chain - name of the chain in the PDB model to examine

    trimer - chain index of middle residue in trimer (NOT CHILD_LIST INDEX)

    score_table - the particular table to lookup against

    """    
    from pro_angle import calc_dihedral
    from math import floor
    
    try:
        (phi,psi) = calc_dihedral(chain, res_id) 
        indx = int(floor(phi/10)+18)
        indy = int(floor(psi/10)+18)

        score = float(score_table[indy][indx])                
    except ValueError:            
        #print "ValueError: asked for score of non-scorable residue"
        score = 0
    return score

def load_scores(score_file):
    """
    Helper function, loads score tables from file

    arguments:

    score_file: location of score table to load

    Returns score table in python format    
    """
    from copy import copy
    from string import atof
    try:
        table = open(score_file)
    except IOError:
        return 0
    tableline = table.readline().split()    
    probx = [0 for i in xrange(36)] #this will be x index
    proby = [0 for i in xrange(36)] #this will be y index    

    for row_counter in range(36):
        for column_counter in range(36):
            probx[column_counter] = atof(tableline[column_counter])
        tableline = table.readline().split()
        proby[row_counter] = copy(probx)
    table = copy(proby)

    return table
        
#main

import sys
import string

if len(sys.argv) > 1:
    PDBfile = sys.argv[1]
else:
    PDBfile = 'error'
if len(sys.argv) > 2:
    start = string.atoi(sys.argv[2])
else:
    start = 1
if len(sys.argv) > 3:
    end = string.atoi(sys.argv[3])    
else:
    end = 9999
if len(sys.argv) > 4:
    smooth = sys.argv[4]
else:
    smooth = 'n'    ###DEFAULT NOT SMOOTHED
if len(sys.argv) > 5:
    mutated = sys.argv[5]
else:
    mutated = ' '

a=score(PDBfile, start, end, smooth, mutated)
print a

