from glob import glob
from sys import argv
import sqlite3 as sql

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
    try:
        with sql.connect(db_file) as conn:
            cur = conn.cursor()
            cur.execute("DROP TABLE IF EXISTS AdjacencyGraph;")
            cur.execute("CREATE TABLE AdjacencyGraph "
                    "(Prefix TEXT PRIMARY KEY, Suffixes TEXT);")
            count = 0
            redundancies = 0
            for file_name in file_names:
                with open(file_name) as f:
                    for line in f:
                        count += 1
                        if count%pattern[0] != pattern[1]:
                            continue
                        line = line.strip()
                        lines = (line[:-1], line[1:])
                        cur.execute("SELECT Suffixes FROM AdjacencyGraph "
                                "WHERE Prefix=?", (lines[0],))
                        data = cur.fetchone()
                        if data:
                            if lines[1] not in data[0]:
                                cur.execute("UPDATE AdjacencyGraph SET "
                                        "Suffixes=?||','||?  WHERE Prefix=?",
                                        (data[0], lines[1], lines[0]))
                                redundancies += 1
                        else:
                            uniques += 1
                            cur.execute("INSERT INTO AdjacencyGraph "
                                    "VALUES(?,?);", lines)
            f = open(status_file, 'a')
            f.write("Initial De Bruijn database populated\n")
            f.write("# of instances of same prefix mapped to distinct suffixes: ")
            f.write(str(redundancies)+"\n")
            f.close()
    finally:
        conn.close()
    return uniques

def traverse_graph(db_file, entries, contig_file, min_length,
            status_file, report_frequency=1000000):
    f = open(contig_file, 'w')
    f.write("")
    f.close()
    try:
        with sql.connect(db_file) as conn:
            c1 = conn.cursor()
            c2 = conn.cursor()
            c1.execute("SELECT * FROM AdjacencyGraph")
            count = 0
            for row in c1:
                count += 1
                if count%report_frequency == 0:
                    f = open(status_file, 'a')
                    f.write(str(100*count/entries)+"% complete\n")
                    f.close()
                if ',' in row[1]:
                    f = open(status_file, 'a')
                    f.write("Ambiguous pairing ignored at entry "+str(count)+
                            ", sequence "+row[0]+"\n")
                    f.close()
                else:
                    prefixes = set(row)
                    prefix = row[1]
                    contig = row[0][0:2]
                    while True:
                        c2.execute("SELECT Suffixes FROM AdjacencyGraph "
                                "WHERE Prefix=?", (prefix,))
                        data = c2.fetchone()
                        if not data:
                            break
                        elif ',' in data[0]:
                            f = open(status_file, 'a')
                            f.write("Sequence mapped to ambiguously " +
                                    "paired prefix "+prefix+". Aborted.\n")
                            f.close()
                            break
                        prefix = data[0]
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
    finally:
        conn.close()

if len(argv) >= 2:
    lines = []
    files = glob(argv[1]) if len(argv) == 2 else argv[1:]
    entries = construct_graph(files, "adjacency_graph.db", "status.txt", (4,2))
    traverse_graph("adjacency_graph.db", entries, "contigs.txt", 160, "status.txt", 500000)

