sample_bam=left_5.bam
reference_fasta=reference/hg19.with.mt.fasta
design_bed=Twist_DNA_ST/pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.bed
output_vcf=

samtools index $sample_bam

gatk --java-options "-Xmx4g" Mutect2 --input $sample_bam --output $output_vcf --reference $reference_fasta --intervals $design_bed
