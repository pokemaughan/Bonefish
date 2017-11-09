from Bio import SeqIO
count = 0
good_reads = (rec for rec in \
              SeqIO.parse("pac-bf_lane3_reads-ACTAGGAG.fastq", "fastq")\
              if min(rec.letter_annotations["phred_quality"])>=20)
count = SeqIO.write(good_reads, "good_quality.fastq","fastq")
print("Saved %i reads" % count)
