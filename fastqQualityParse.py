from Bio import SeqIO

def qualityParse(filename,quality):
    count = 0
    good_reads = (rec for rec in \
              SeqIO.parse(filename, "fastq")\
              if min(rec.letter_annotations["phred_quality"])>=quality)
    count = SeqIO.write(good_reads, "good_quality.fastq","fastq")
    print("Saved %i reads" % count)

def main():
    qualityParse("pac-bf_lane3_reads-AGATACAA.fastq",30)

if __name__=='__main__':
    main()



