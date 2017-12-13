import difflib
import time

start_time = time.time()

def getOverlaps(file):
    superread = ""
    superreads = list()
    matches = 0
    matcher = difflib.SequenceMatcher(None)
    for line in range(len(file)):
        matcher.set_seqs(file[line],file[line])
        biggestmatch = matcher.find_longest_match(0,len(file[line]),0,len(file[line]))
        read1 = line
        read2 = line
        for line2 in range(line,len(file)):
            matcher.set_seqs(file[line],file[line2])
            match = matcher.find_longest_match(0,len(file[line]),0,len(file[line2]))
            if(match.size > 50 and match.size < 150 and match.b==0):
                if(match.a > biggestmatch.a):
                    biggestmatch = match
                    read2 = line2
                    
                
        print(biggestmatch)
        print("%d %d" %(read1,read2))
        if (read1 != read2):
##            print(file[read1])
##            s = "";
##            for x in range(biggestmatch.a):
##                s +=" "
##            matching = file[read2]
##            strings = "%s%s" %(s,matching)
##            print(strings)
##            print ("-----------------------------------------------------------")
            part1 = (file[read1][:biggestmatch.a])
            #print (part1)
            part2 = (file[read2])
            #print (part2)
            superread = "%s%s" %(part1.strip(),part2.strip())
            #print (file[read1])
            print (superread)
            superreads.append(superread)
        
        
    print(matches)        
    return superreads

def main():
    count=0;
    sums =0;
    contigs = list();
    with open('new_contigs.txt') as input_data:
        for line in input_data:
            count+=1

            #strings = "%d\t%d" %(count,len(line))
            #print(strings)
            contigs.append(line)

        contigs.sort(key = len)
        medianpos = int(len(contigs)/2)
        N50 = len(contigs[medianpos])
        longest = len(max(contigs,key=len))

    overlaps = getOverlaps(contigs)

    count2 = len(overlaps)
    overlaps.sort(key = len)
    medianpos = int(len(overlaps)/2)
    N502 = len(overlaps[medianpos])
    longest2 = len(max(overlaps,key=len))

    # Print and save the answer.
    print(count)
    print (N50)
    print (longest)
    with open('ContigAnalysis.txt', 'w') as output_data:
        count_s = "Count:\t %d %d\n" %(count,count2)
        N50_s = "N50:\t %d %d\n" %(N50,N502)
        longest_s = "Longest: %d %d\n" %(longest,longest2)
        output_data.write(str(count_s))#(", ".join(map(str,edges)))
        output_data.write(str(N50_s))
        output_data.write(str(longest_s))
        output_data.write(("--- %s seconds ---" % (time.time() - start_time)))

if __name__=='__main__':
    main()
