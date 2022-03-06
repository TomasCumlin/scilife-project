import sys
import pysam

threshold=sys.argv[2]

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

counter = 0
textfile = open("filtered_names.txt","w")

for read in samfile:
	if "MD" in str(read.tags[1]):
		counter += str(read.tags[1]).count("C")
		counter += str(read.tags[1]).count("c")
		counter += str(read.tags[1]).count("G")
		counter += str(read.tags[1]).count("g")
		counter += str(read.tags[1]).count("A")
		counter += str(read.tags[1]).count("a")
		counter += str(read.tags[1]).count("T")
		counter += str(read.tags[1]).count("t")
		if counter > int(threshold):
			textfile.write(str(read.qname) + "\n")
		counter = 0

textfile.close()
samfile.close()
