#take mscores and make them plottable

big = open('mscores')
b = big.readline()
new = open('mscores_form', 'w')
while b:
    num = b[42:51]
    if num != '0\n' and num != '0' and num[0] != '0':
        new.write(num)
        new.write('\n')
    b = big.readline()
new.close()
        
