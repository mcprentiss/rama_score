#parse through all of mike's data, split it into residue specific data
#/home/mpretiss/proteins/dist/data
import os
import commands


aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
       'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
       'ARG', 'HIS']

os.chdir('/home/mprentis/proteins/dist/data')

list = commands.getoutput('ls *.dat').split()


for aa in aas:
    print aa
    os.chdir('/home/marciovm/proteins/')
    file = open(aa + '.dat', 'w')
    os.chdir('/home/mprentis/proteins/dist/data')
    for pdb in list:
            readin = open(pdb, 'r')
            count = 0

            while count < 20:
                line = readin.readline()                            
                if len(line) < 2:                
                    pass
                else:
                    code = line[9]+line[10]+line[11]
                    if code == aa:
                        match = 0

                        line = readin.readline()
                        while len(line) > 2:
                           
                            file.write(line)
                            line = readin.readline()                
                    else:
                        while len(line) > 2:
                            line = readin.readline()
                count += 1
    file.close()
                    
                    
        
    



