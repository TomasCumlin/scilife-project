# Quality report (and trimming if set to that)

fastp -i $R1 -I $R2 -t $tail $front -o $batch.R1.fq.gz -O $batch.R2.fq.gz

bwa mem -t 4 -R '@RG\tID:$sample\tSM:$sample' -v 1 $reference $batch.R1.fq.gz $batch.R2.fq.gz | samtools sort -@ 4 -o $batch.bam

samtools index $batch.bam

sh markduplicates_code.sh

sh hsmetrics_duplicatemetrics_code.sh

$multiqc_command
