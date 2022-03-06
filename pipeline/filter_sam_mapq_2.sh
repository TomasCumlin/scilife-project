# python script filter_sam_mapq.py  must be located in the same folder as this sh-script.

bamfile=marked_duplicates_$batch.bam
outfile=mapq_${threshold_mapq}_$batch.bam

# get text-file with read names from  bamfile
python filter_sam_mapq_2.py $bamfile $threshold_mapq


picard FilterSamReads \
    INPUT=$bamfile \
    OUTPUT=$outfile \
    READ_LIST_FILE=mapq_low_list.txt \
      FILTER=excludeReadList


outfile2=mapqout_${threshold_mapq}_$batch.bam


picard FilterSamReads \
    INPUT=$bamfile \
    OUTPUT=$outfile2 \
    READ_LIST_FILE=mapq_low_list.txt \
      FILTER=includeReadList


#rm mapq_low_list.txt
