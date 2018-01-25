### Plot trimming statistics
# Simple R script to visualise the effect of trimming on the FastQ files

### Read data in

tbl <- read.table("./trim/trimOutput/outputTrimStats.txt")


### Calculate percentage trimming

PercentTrimmedReads <- tbl[,3]/tbl[,1]*100
PercentQualityTrimmed <- tbl[,4]/tbl[,2]*100
PercentTrimmedBases <- tbl[,5]/tbl[,2]*100


### Bind statistics to dataframe, apply averaging

stats <- data.frame(PercentTrimmedReads, PercentQualityTrimmed,PercentTrimmedBases)
avr <- apply(stats,2,mean)


### Generate PDF plots

pdf("./trim/trimOutput/TrimStats_PercTrimmedReads.pdf",8,8)
  barplot(stats[,1], xlab="FastQ files", ylab="% Trimmed Reads")
dev.off()

pdf("./trim/trimOutput/TrimStats_PercQualityTrimmed.pdf",8,8)
  barplot(stats[,2], xlab="FastQ files", ylab="% Quality Trimmed")
dev.off()

pdf("./trim/trimOutput/TrimStats_PercTrimmedBases.pdf", 8,8)
  barplot(stats[,3], xlab="FastQ files", ylab="% Trimmed Bases")
dev.off()

pdf("./trim/trimOutput/TrimStats_Averages.pdf",8,8)
  barplot(avr, ylab="Percentage", main="Percentage trimming averaged across all files")
dev.off()
