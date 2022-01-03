#bwa index marked_duplicates_filtered.bam

sample_bam=marked_duplicates.bam
reference_fasta=reference/hg19.with.mt.fasta
design_bed=Twist_DNA_ST/pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.bed

gatk --java-options "-Xmx4g" Mutect2 --input $sample_bam --output unfiltered_bam.vcf --reference $reference_fasta --intervals $design_bed
