samtools index $filtered_bam

gatk --java-options "-Xmx4g" Mutect2 --input $filtered_bam --output $batch.vcf --reference $reference --intervals $intervals

python vc_check_updated.py $batch.vcf HD832_variants.txt $batch
