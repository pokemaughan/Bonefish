output = open("reduce_status.txt", 'w')
output.write("")
output.close()

with open("contigs.txt") as f:
    lines = [line.strip() for line in f.readlines()]

def substringSieve(string_list):
    output = open("reduce_status.txt", 'w')
    output.write("Began substring sieve\n")
    output.close()
    string_list.sort(key=lambda s: len(s), reverse=True)
    out = []
    count = 0
    length = len(string_list)
    for s in string_list:
        if count%100 == 0:
            output = open("reduce_status.txt", 'w')
            output.write(str(count/length) + "% parsed\n")
            output.close()
        if not any([s in o for o in out]):
            out.append(s)
        count += 1
    return out

f = open("reduced_contigs.txt", 'w')
reduced_lines = substringSieve(lines)
count = 0
length = len(reduced_lines)
for line in reduced_lines:
    if count%100 == 0:
        output = open("reduce_status.txt", 'w')
        output.write(str(count/length) + "% reported\n")
        output.close()
    f.write(line + "\n")
    count += 1
f.close()