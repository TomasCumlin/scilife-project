# multiqc duplicates

dup_table = read.table("multiqc_picard_dups.txt", header = TRUE, sep = "", dec = ".")

# hsmetrics

hs_table = read.table("multiqc_picard_HsMetrics.txt", header = TRUE, sep = "", dec = ".", fill = TRUE)
