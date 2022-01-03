picard FilterSamReads \
    INPUT=marked_duplicates.bam \
    OUTPUT=marked_duplicates_filtered.bam \
    READ_LIST_FILE=text_file.txt \
      FILTER=excludeReadList
