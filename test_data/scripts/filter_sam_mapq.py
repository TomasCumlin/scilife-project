import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

textfile = open("mapq_low_list.txt","w")

for i in samfile:
	if i.mapping_quality < 30:
		textfile.write(str(i.qname) + "\n")

textfile.close()
samfile.close()
