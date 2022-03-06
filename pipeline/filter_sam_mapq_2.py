import sys
import pysam

threshold=int(sys.argv[2])

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

mapq_textfile = open("mapq_low_list.txt","w")

for i in samfile:
	if i.mapping_quality <= threshold:
		mapq_textfile.write(str(i.qname) + "\n")

mapq_textfile.close()
samfile.close()
