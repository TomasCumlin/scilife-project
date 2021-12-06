import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

counter = 0
textfile = open("text_file.txt","w")

for read in samfile:
	if "MD" in str(read.tags[1]):
		counter += str(read.tags[1]).count("C")
		counter += str(read.tags[1]).count("G")
		counter += str(read.tags[1]).count("A")
		counter += str(read.tags[1]).count("T")
		if counter > 20:
			textfile.write(str(read.qname) + "\n")
		counter = 0

textfile.close()
samfile.close()
