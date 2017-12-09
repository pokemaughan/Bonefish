from glob import glob
from sys import argv
import sqlite3 as sql

all_reads = []

def construct_graph(file_names, status_file):
    """Constructs de Bruijn graph of reads, formatted as relational database.

    Parameters:
        file_names: list of string names of .fastq files with reads to process.
        db_file: string name of the database file to be written to
        status_file: string name of file to write status updates in

    Returns: number of lines (or unique prefixes) in the database"""
    f = open(status_file, 'w')
    f.write("")
    f.close()
    uniques = 0
    adjacency_graph = {}
    first_to_last = {}
    last_to_first = {}
    count = 0
    redundancies = 0
    current_first = ""
    identifier = ""
    for file_name in file_names:
        with open(file_name) as f:
            paired = False
            for line in f:
                count += 1
                if count%4 not in [1,2]:
                    continue
                elif count%4 == 1:
                    if line[:-2] == identifier[:-2]:
                        paired = True
                    else:
                        identifier = line
                    continue
                all_reads.append(line)
                line = line.strip()
                pairing = (line[:-1], line[1:])
                if paired:
                    if current_first not in first_to_last:
                        first_to_last[current_first] = set([pairing[0]])
                    else:
                        first_to_last[current_first].add(pairing[0])
                    if pairing[0] not in last_to_first:
                        last_to_first[pairing[0]] = set([current_first])
                    else:
                        last_to_first[pairing[0]].add(current_first)
                else:
                    current_first = pairing[0]
                paired = False
                if pairing[0] in adjacency_graph:
                    if pairing[1] not in adjacency_graph[pairing[0]]:
                        adjacency_graph[pairing[0]] += "," + pairing[1]
                        redundancies += 1
                else:
                    uniques += 1
                    adjacency_graph[pairing[0]] = pairing[1]
    f = open(status_file, 'a')
    f.write("Initial De Bruijn database populated\n")
    f.write("# of instances of same prefix mapped to distinct suffixes: ")
    f.write(str(redundancies)+"\n")
    f.close()
    return uniques, adjacency_graph, first_to_last, last_to_first

def traverse_graph(adjacency_graph, forward_pairs, backward_pairs, entries, contig_file, min_length,
            status_file, report_frequency=1000000):
    f = open(contig_file, 'w')
    f.write("")
    f.close()
    count = 0
    for prefix in adjacency_graph:
        count += 1
        if count%report_frequency == 0:
            f = open(status_file, 'a')
            f.write(str(100*count/entries)+"% complete\n")
            f.close()
        if ',' in adjacency_graph[prefix]:
            print("\n*****")
            print("\tPrefix:", prefix)
            print("\tSuffixes:")
            for suffix in str(adjacency_graph[prefix]).split(','):
                print('\t'+suffix)
                pairs = []
                if suffix in forward_pairs:
                    pairs.extend(forward_pairs[suffix])
                if suffix in backward_pairs:
                    pairs.extend(backward_pairs[suffix])
                for pair in pairs:
                    print('\t\t', pair, sep='')
            pairs = []
            if prefix in forward_pairs:
                pairs.extend(forward_pairs[prefix])
            if prefix in backward_pairs:
                pairs.extend(backward_pairs[prefix])
            print("\tPairs:")
            for pair in pairs:
                print(pair)
            print("*****\n")
            f = open(status_file, 'a')
            f.write("Ambiguous pairing ignored at entry "+str(count)+
                    ", sequence "+prefix+"\n")
            f.close()
        else:
            prefixes = set(adjacency_graph[prefix])
            prefixes.add(prefix)
            contig = prefix[0:2]
            prefix = adjacency_graph[prefix]
            while True:
                if prefix not in adjacency_graph:
                    break
                elif ',' in adjacency_graph[prefix]:
                    f = open(status_file, 'a')
                    f.write("Sequence mapped to ambiguously " +
                            "paired prefix "+prefix+". Aborted.\n")
                    f.close()
                    break
                prefix = adjacency_graph[prefix]
                if prefix in prefixes:
                    f = open(status_file, 'a')
                    f.write("Infinite loop detected. Path ignored.\n")
                    f.close()
                    break
                prefixes.add(prefix)
                contig += prefix[0]
            contig += prefix[1:]
            if len(contig) >= min_length:
                f = open(contig_file, 'a')
                f.write(contig + "\n")
                f.close()

if len(argv) >= 2:
    lines = []
    files = glob(argv[1]) if len(argv) == 2 else argv[1:]
    entries, adjacency_graph, forward_pairs, backward_pairs = construct_graph(files, "status.txt")
    traverse_graph(adjacency_graph, forward_pairs, backward_pairs, entries, "contigs.txt", 152, "status.txt", 500000)

