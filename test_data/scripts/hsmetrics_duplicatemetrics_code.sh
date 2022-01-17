input_file = 

java -jar picard.jar CollectHsMetrics COVERAGE_CAP=5000 I=$input_file O=hs_metrics_left_5.txt R=/home/tomas/data/reference/hg19.with.mt.fasta BAIT_INTERVALS=pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.interval_list  TARGET_INTERVALS=pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.interval_list

java -jar picard.jar CollectDuplicateMetrics INPUT=$input_file M=collectedduplicate_left_5.metrics
