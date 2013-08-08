#take mscores3 and make them plottable

big = open('mscores3_1')
new = open('mscores3_form', 'w')

big2 = open('mscores3_2')
big3 = open('mscores3_3')

list = [big, big2, big3]

for set in list:
    b = set.readline()
    while b:
        start = 9
        bump = 0
        if b[9] == ' ':
            bump = 1       
        num = b[(9 + bump):(15 + bump)]
        if num != '0\n' and num != '0' and num[0] != '0':
            new.write(num)
            new.write('\n')
        b = set.readline()

new.close()
        
