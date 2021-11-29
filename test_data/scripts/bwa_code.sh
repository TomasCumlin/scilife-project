sample=sample

bwa mem -t 4 -R '@RG\tID:${sample}\tSM:${sample}' -v 1 /home/tomas/data/reference/hg19.with.mt.fasta /home/tomas/data/Twist_DNA_ST/VAL-01_R1.fastq.gz /home/tomas/data/Twist_DNA_ST/VAL-01_R2.fastq.gz | samtools sort -@ 4 -o ${sample}.bam
