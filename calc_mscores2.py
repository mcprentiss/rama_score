#! /usr/bin/env python

import os
import commands
import mscore3
import sys
import string

filename = sys.argv[1]

start = sys.argv[2]
end = sys.argv[3]

file_out = open(filename, 'a')

os.chdir('/home/marciovm/proteins/bigdist/')
cmd  = 'ls *.*'
file_list = commands.getoutput(cmd).split()

os.chdir('/home/marciovm/proteins')
for index in range(string.atoi(start), string.atoi(end)):
    file = file_list[index]
    print "calculating score of " + file + ' ' + str(index)
    scr = mscore3.score('/home/marciovm/proteins/bigdist/' + file)  
    file_out.write(file + ' ' + str(scr) + '\n')    
file_out.close()
    
    
