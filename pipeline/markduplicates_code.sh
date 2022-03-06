java -jar picard.jar MarkDuplicates \
      I=$batch.bam \
      O=marked_duplicates_$batch.bam \
      M=marked_duplicates_metrics_$batch.txt

samtools index marked_duplicates_$batch.bam
