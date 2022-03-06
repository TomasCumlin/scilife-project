python sam_test.py marked_duplicates_$batch.bam $threshold

input_bam=marked_duplicates_$batch.bam
output_bam=filteredout_${threshold}_$batch.bam

picard FilterSamReads \
    INPUT=$input_bam \
    OUTPUT=$output_bam \
    READ_LIST_FILE=filtered_names.txt \
      FILTER=includeReadList

output_bam2=filtered_${threshold}_$batch.bam

picard FilterSamReads \
    INPUT=$input_bam \
    OUTPUT=$output_bam2 \
    READ_LIST_FILE=filtered_names.txt \
      FILTER=excludeReadList

