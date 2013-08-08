#script to get representative proline angles from protlist_alpha
import os
import commands
from pro_angle import pro_angle

file = open('protlist_alpha', 'r')
list = file.read().split()
file_out = open('PROLINES', 'w')

for pdb in list:
    print "processing... " + str(pdb)
    os.chdir('/rsyncpdb/data/structures/all/pdb/')
    fourletter =  pdb[0]+pdb[1]+pdb[2]+pdb[3]
    filename = 'pdb' + fourletter.lower() + '.ent.Z'
    parts = filename.split('.')
    cmd = 'cp ' + filename + ' ~/proteins/dist'
    print cmd
    commands.getoutput(cmd)
    os.chdir('/home/marciovm/proteins/dist/')
    cmd2= 'mv ' + filename + ' ' + parts[0] + pdb[4] + '.ent.Z'
    print cmd2
    commands.getoutput(cmd2)
    commands.getoutput('rm ' + filename)

