# python script filter_sam_mapq.py  must be located in the same folder as this sh-script.

bamfile=sample.bam
outfile=filtered.bam

python filter_sam_mapq.py $bamfile


picard FilterSamReads \
    INPUT=$bamfile \
    OUTPUT=$outfile \
    READ_LIST_FILE=mapq_low_list.txt \
      FILTER=excludeReadList

rm mapq_low_list.txt
