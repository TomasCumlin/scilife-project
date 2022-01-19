input=marked_duplicates.bam
reference=reference/hg19.with.mt.fasta

java -Xmx4g -jar /home/tomas/miniconda3/pkgs/fgbio-1.5.0-hdfd78af_0/share/fgbio/fgbio.jar ErrorRateByReadPosition -i $input -r $reference -m 0 --collapse=False

