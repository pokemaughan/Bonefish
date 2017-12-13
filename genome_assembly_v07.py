from glob import glob
from sys import argv
import time

def construct_graph(file_names, status_file, pattern=(1,1)):
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
    redundancies = 0
    for file_name in file_names:
        count = 0
        with open(file_name) as f:
            for line in f:
                count += 1
                if count%pattern[0] != pattern[1]:
                    continue
                line = line.strip()
                prefix, suffix = line[:-1], line[1:]
                if prefix in adjacency_graph:
                    if suffix not in adjacency_graph[prefix]:
                        adjacency_graph[prefix] += "," + suffix + "|1"
                        redundancies += 1
                    else:
                        suffixes = adjacency_graph[prefix].split(',')
                        for i in range(0, len(suffixes)):
                            if suffix in suffixes[i]:
                                data = suffixes[i].split('|')
                                string, number = data[0], int(data[1])
                                suffixes[i] = string + "|" + str(number+1)
                                break
                        adjacency_graph[prefix] = ','.join(suffixes)
                else:
                    uniques += 1
                    adjacency_graph[prefix] = suffix+"|1"
        f = open(status_file, 'a')
        f.write("Finished reading " + file_name + "\n")
        f.close()
    f = open(status_file, 'a')
    f.write("Initial De Bruijn database populated\n")
    f.write("# of instances of same prefix mapped to distinct suffixes: ")
    f.write(str(redundancies)+"\n")
    f.close()
    return uniques, adjacency_graph

def traverse_graph(adjacency_graph, entries, contig_file, min_length,
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
        if "," in adjacency_graph[prefix]:
            pass
            # f = open(status_file, 'a')
            # f.write("Ambiguous pairing ignored at entry "+str(count)+
            #         ", sequence "+prefix+"\n")
            # f.close()
        else:
            suffix = list(adjacency_graph[prefix].split('|'))[0]
            prefixes = set([suffix])
            prefixes.add(prefix)
            contig = prefix[0:2]
            prefix = suffix
            while True:
                suffix = ""
                if prefix not in adjacency_graph:
                    break
                elif "," in adjacency_graph[prefix]:
                    f = open(status_file, 'a')
                    f.write("Counts for ambiguous pairing: ")
                    best_suffix = ""
                    best_count = -1
                    suffixes = adjacency_graph[prefix].split(',')
                    for suffix in suffixes:
                        data = suffix.split('|')
                        this_suffix, this_count = data[0], int(data[1])
                        f.write(str(this_count)+" ")
                        if this_count > best_count:
                            best_suffix = this_suffix
                            best_count = this_count
                    f.write('\n')
                    f.close()
                    suffix = best_suffix
                else:
                    suffix = list(adjacency_graph[prefix].split('|'))[0]
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
    entries, adjacency_graph = construct_graph(files, "status.txt", (4,2))
    traverse_graph(adjacency_graph, entries, "contigs.txt", 152, "status.txt", 500000)
    print(time.time() - start_time)

