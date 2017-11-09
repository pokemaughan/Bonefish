def string_from_paired_composition(k, d, paired_reads):
    mapping = {}
    for p in paired_reads:
        for r in paired_reads:
            if (p[0][:-1] == r[0][1:] and p[1][:-1] == r[1][1:]):
                mapping[(p[0],p[1])] = (r[0], r[1])
    missing = [x for x in paired_reads if (x[0], x[1]) not in mapping][0]
    mapping = {}
    for p in paired_reads:
        for r in paired_reads:
            if (p[0][1:] == r[0][:-1] and p[1][1:] == r[1][:-1]):
                mapping[(p[0],p[1])] = (r[0], r[1])
    front = missing[0]
    back = missing[1]
    node = missing
    while node in mapping:
        front += mapping[node][0][-1]
        back += mapping[node][1][-1]
        node = mapping[node]
    string = "_"*(2*k+d+len(paired_reads)-1)
    string = front + string[-1*(len(string)-k-len(paired_reads)+1-1):]
    string = string[:len(string)-k-len(paired_reads)+2] + back
    return string

with open('rosalind_ba3j.txt') as f:
    lines = [line.strip() for line in f.readlines()]
k, d = lines[0].split(' ')
f = open('output.txt', 'w')
tuples=[]
for L in [x.split('|') for x in lines[1:]]:
    tuples.append((L[0],L[1]))
f.write(string_from_paired_composition(int(k), int(d), tuples))
