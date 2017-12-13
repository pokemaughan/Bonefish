from glob import glob
from sys import argv
import time

def construct_graph(file_names, db_file, status_file, pattern=(1,1)):
    """Constructs de Bruijn graph of reads, formatted as relational database.

    Parameters:
        file_names: list of string names of files with reads to process.
        db_file: string name of the database file to be written to
        status_file: string name of file to write status updates in
        pattern: 2-tuple of ints (n, k). For every set of n lines in the input files,
                all but the kth line will be ignored. Default (1,1) reads all lines.

    Returns: number of lines (or unique prefixes) in the database"""
    f = open(status_file, 'w')
    f.write("")
    f.close()
    uniques = 0
    adjacency_graph = {}
    count = 0
    redundancies = 0
    for file_name in file_names:
        with open(file_name) as f:
            for line in f:
                count += 1
                if count%pattern[0] != pattern[1]:
                    continue
                line = line.strip()
                prefix, suffix = line[:-1], line[1:]
                if prefix in adjacency_graph:
                    if suffix not in adjacency_graph[prefix]:
                        adjacency_graph[prefix][suffix] = 1
                        redundancies += 1
                    else:
                        adjacency_graph[prefix][suffix] += 1 #"that specific suffix is in the graph and is being hit by the search again"
                else:
                    uniques += 1
                    adjacency_graph[prefix] = {suffix:1}
        f = open(status_file, 'a')
        f.write("Finished reading " + file_name + "\n")
        f.close()
    f = open(status_file, 'a')
    f.write("Initial De Bruijn database populated\n")
    f.write("# of instances of same prefix mapped to distinct suffixes: ")
    f.write(str(redundancies)+"\n")
    f.close()
    return uniques, adjacency_graph

def traverse_graph(adjacency_graph, db_file, entries, contig_file, min_length,
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
        if len(adjacency_graph[prefix]) > 1:
            pass
            # f = open(status_file, 'a')
            # f.write("Ambiguous pairing ignored at entry "+str(count)+
            #         ", sequence "+prefix+"\n")
            # f.close()
        else:
            prefixes = set(adjacency_graph[prefix].keys())
            prefixes.add(prefix)
            contig = prefix[0:2]
            prefix = list(adjacency_graph[prefix].keys())[0]
            while True:
                suffix = ""
                if prefix not in adjacency_graph:
                    break
                elif len(adjacency_graph[prefix]) > 1:
                    f = open(status_file, 'a')
                    f.write("Counts for ambiguous pairing: ")
                    best_suffix = ""
                    best_count = -1
                    for suffix in adjacency_graph[prefix]:
                        this_count = adjacency_graph[prefix][suffix].count
                        f.write(str(this_count)+" ")
                        if this_count > best_count:
                            best_suffix = suffix
                            best_count = this_count
                    f.write('\n')
                    f.close()
                    suffix = best_suffix
                else:
                    suffix = list(adjacency_graph[prefix].keys())[0]
                prefix = suffix
                if prefix in prefixes:
                    # f = open(status_file, 'a')
                    # f.write("Infinite loop detected. Path ignored.\n")
                    # f.close()
                    break
                prefixes.add(prefix)
                contig += prefix[0]
            contig += prefix[1:]
            if len(contig) >= min_length:
                f = open(contig_file, 'a')
                f.write(contig + "\n")
                f.close()

if len(argv) >= 2:
    start_time = time.time()
    lines = []
    files = glob(argv[1]) if len(argv) == 2 else argv[1:]
    entries, adjacency_graph = construct_graph(files, "adjacency_graph.db", "status.txt", (4,2))
    traverse_graph(adjacency_graph, "adjacency_graph.db", entries, "contigs.txt", 175, "status.txt", 500000)
    print (time.time() - start_time)

