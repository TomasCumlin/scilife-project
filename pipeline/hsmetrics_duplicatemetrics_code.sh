input_file=marked_duplicates_$batch.bam
output_txt=hs_metrics_$batch.txt
metrics_out=collectedduplicate_$batch.metrics

java -jar picard.jar CollectHsMetrics COVERAGE_CAP=5000 I=$input_file O=$output_txt R=$reference BAIT_INTERVALS=$intervals  TARGET_INTERVALS=$intervals

java -jar picard.jar CollectDuplicateMetrics INPUT=$input_file M=$metrics_out

