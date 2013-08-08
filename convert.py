file = open('PROLINES', 'r')
out = open('PRO', 'w')

a = file.readline()
while a:
    if a[0] != '#':
        out.write(a)
    a = file.readline()

out.close()
        
        
