from Bio import SeqIO
count = 0
for rec in SeqIO.parse("pac-bf_lane3_reads-ACTAGGAG.fastq", "fastq"):
    count += 1
print("%i reads" % count)
