library(tools)

# Import output
tbl <- read.table("./StarBAM/Log/outputMappingStats.txt", header=T, row.names=1)
rownames(tbl) <- basename(file_path_sans_ext(file_path_sans_ext(rownames(tbl))))

# Read lengths
readLength <- tbl[,c(2,4)]
pdf("./StarBAM/Log/readLength.pdf",8,8)
  boxplot(readLength, main = "Read length")
dev.off()

# Mapping rate

readMappingPerc <- data.frame(tbl[,3]/tbl[,1]*100, tbl[,6]/tbl[,1]*100, tbl[,8],tbl[,9],tbl[,10])
colnames(readMappingPerc) <- c("% Uniquely mapped", 
                               "% Multimapped", 
                               "% Unmapped - mismatches",
                               "% Unmapped - short read",
                               "% Unmapped - other")

pdf("./StarBAM/Log/mappingRates.pdf",12,8)
  par(mar=c(11,4,4,2))
  boxplot(readMappingPerc, las=2, main="Boxplot of STAR mapping rates", ylab="Percentage")
dev.off()

# Raw counts
readRaw <- tbl[,c(1,3,6)]

pdf("./StarBAM/Log/rawReadCounts.pdf", 12,8)
  par(mar=c(9,4,6,2))
  barplot(t(readRaw), 
          las=2, 
          main="Barplot of total read counts", 
          beside=T, 
          col=c("Black","Grey", "White"),
          ylim=c(0,max(readRaw)*1.5))
  legend("topright",
          cex=0.8,
          legend=c("Total reads", "Uniquely mapped reads", "Multimapped reads"), 
          fill=c("Black","Grey", "White")
)
