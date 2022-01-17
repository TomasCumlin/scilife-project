inFile=marked_duplicates.bam
outFile=left_5.bam
left=5
right=0

/home/tomas/miniconda3/bin/bam trimBam $inFile $outFile -L $left -R $right
