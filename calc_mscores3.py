#! /usr/bin/env python

import os
import commands
import mscore3py
import sys
import string

filename = sys.argv[1]

start = sys.argv[2]
end = sys.argv[3]

file_out = open(filename, 'w')

os.chdir('/home/marciovm/proteins/bigdist/')
cmd  = 'ls *.*'
file_list = commands.getoutput(cmd).split()

os.chdir('/home/marciovm/proteins')
for index in range(string.atoi(start), string.atoi(end)):
    file = file_list[index]
    print "calculating score of " + file + ' ' + str(index)
    scr = mscore3py.score('/home/marciovm/proteins/bigdist/' + file)  
    file_out.write(file + ' ' + str(scr) + '\n')    
file_out.close()
    
    
