#compile proline phi chi psi angles
import os
import commands
from pro_angle import pro_angle

os.chdir('/home/marciovm/proteins/dist')
list = commands.getoutput('ls').split()

file_out = open('PROLINES', 'w')
count = 0

for pdb in list:
    commands.getoutput("gunzip " + pdb)
    if pdb[7] == '_':
        angle_list=pro_angle(pdb.rstrip('.Z'), 'PRO', ' ', 0, 0, 'y')
    else:
        angle_list=pro_angle(pdb.rstrip('.Z'), 'PRO', pdb[7], 0, 0, 'y')

    print 'writing output... '
    
    file_out.write('#pdb: ' +  pdb + '\n')
    
    for tuple in angle_list:
        if tuple:
            file_out.write(str(tuple[0]).zfill(6) + '  '
                           +str(tuple[1]).zfill(6)+ '  '
                           +str(tuple[2]).zfill(6)+ '  '
                           +str(tuple[3]).zfill(6)+'\n')
            count += 1
        else:
            pass
file_out.write('#occurences of sequence: ' + str(count)+'\n')
file_out.close()
