#bring files listed in input to target directory
#pre define input file and output directory
#need to run in ellison,

import string
import os
import commands

data_dir = '/rsyncpdb/data/structures/all/pdb/'
target_dir ='/home/marciovm/proteins/bigdist/'

file = open('PDB.LIST')
targets = file.read().split()

count = 0
for aa in targets:
    name = string.lower(aa)
    source_filename = data_dir + 'pdb' + name[0:4] + '.ent.Z'
    target_filename = target_dir + name + '.ent.Z'
    cmd = 'cp ' + source_filename + ' ' + target_filename
    commands.getoutput(cmd)
    print 'copied # ' + str(count)
    count += 1
    
