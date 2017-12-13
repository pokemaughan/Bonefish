"""
TODO:
    - Make proper De Bruijn graph
    - Handle and prune bubbles
    - Traverse and find possible contigs
    - Compare to the N50 standard
"""

def string_from_paired_composition(k, d, paired_reads):
    """Assembles paired reads into a single contig
    
    Parameters:
        k (int): length of the kmers at either end of the pair
        d (int): length of unknown DNA segment between paired reads
        paired_reads (list of 2-string tuples): the paired reads
        
    """
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
    if '_' in string:
        raise ValueError("Contig was not assembled correctly. " +
                "Unidentified base pairs exist.")
    return string
