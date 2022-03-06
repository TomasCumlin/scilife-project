# Instructions:
	# Update all variables.
	# Select which codes to use and exlcude by using the #


export batch="test"

export reference=/home/tomas/data/reference/hg19.with.mt.fasta

export intervals=pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.interval_list

export bed=Twist_DNA_ST/pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.bed

# enter fastq-file names
export R1=Twist_DNA_ST/VAL-01_R1.fastq.gz
export R2=Twist_DNA_ST/VAL-01_R2.fastq.gz


# trimming setting (assign 0 if no trimming needed)
export front=0
export tail=10

# MultiQC command.
export multiqc_command="multiqc . --ignore 'Twist_DNA_ST/'"

# <ACTIVATE FASTP>
# <ACTIVATE BWA>
# <ACTIVATE SAMTOOLS>
# <ACTIVATE MULTIQC>
# <ACTIVATE PICARD>
# <ACTIVATE GATK4>
# <ACTIVATE MPILEUP>

# 01_create_bam
# -> Quality check with fastp, trimming if wanted, alignment, mark duplicates (via markduplicates_code.sh), collect metrics and duplicates (via hsmetrics_duplicatemetrics_code.sh), MultiQC

sh 01_create_bam.sh

# Plot Reads Based on # Substitutions present
#python subst_count_2.py marked_duplicates_$batch.bam $batch

# Plot Reads Based on their Mapq-value
#python mapq_counts.py marked_duplicates_$batch.bam $batch

# Filter Reads Based on # substitutions present

# Select # substitutions threshold
export threshold=4

sh filtersamreads.sh

# Filter Reads Based on their Mapq_value

# Select maximum mapq-value
export threshold_mapq=41

#sh filter_sam_mapq_2.sh

# If BAM-file was filtered based on substitutions, use these variables:
export filtered_bam=filtered_${threshold}_$batch.bam
export filteredout_bam=filteredout_${threshold}_$batch.bam

# If BAM-file was filtered based on mapq-value, use these variables:
#export filtered_bam=mapq_${threshold_mapq}_$batch.bam
#export filteredout_bam=mapqout_${threshold_mapq}_$batch.bam

# Create VCF file (followed by a VC-check using vc_check_updated.py)
sh vc_gatk4.sh

# Create Error Profiles

samtools index $filteredout_bam
python pileup_plots_updated.py $batch $filtered_bam $filteredout_bam $reference $bed
