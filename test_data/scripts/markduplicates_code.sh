java -jar picard.jar MarkDuplicates \
      I=sample.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
