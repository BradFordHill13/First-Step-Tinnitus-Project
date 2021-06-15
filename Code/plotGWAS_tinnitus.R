#install.packages("qqman")
#install.packages("BiocManager")
#BiocManager::install("GWASTools")
library(qqman)
library(GWASTools)
setwd("/Users/mtomason/Documents/projects/tinnitus/data_UKBiobank/03_gwas/qqplot")

# INSTRUCTIONS 
# - take GWAS output files from Jura
# OUTPUT
# for each chromo: qqplot and manhattan plot
# genomewide qqplot and manhattan plot

plotPvals <- function(name,pheno,do_qqplot,do_manhattan){
  # rename according to qqman requirements: SNP CHR BP P
  colnames(pheno) <- c("CHR","SNP","BP","P")
  pvalues <- `^`(10,-pheno$P) # transform -log10 p values
  pheno$P <- pvalues
  if(do_qqplot){
    try(GWASTools::qqPlot(pvalues, main=name))
  }
  if(do_manhattan){
    try(qqman::manhattan(pheno, main=name))
  }
}

# aggregate GWAS results from each chromo
gwasResults_allChr__tinnitus <- data.frame()

for (i in c(12:12)){
  write(paste0("processing chromo",i), stdout())
  
  # read input
  gwasResults <- read.table(paste("output_ukb_imp_chr", i,"_v3.txt", sep=""), sep=" ",header=T, stringsAsFactors= F)
  gwasResults <- gwasResults[complete.cases(gwasResults), ] # drop NAs (can happen when maf=1)
  # prepare plot file
  jpeg(file=paste("pvals_chr", i,"_v3.jpg", sep=""), width=6000,height=2000)
  old.par <- par(mfrow=c(1, 2)) # 1 phenos, 2 plots (qqplot/manhatt)
  par(mar=c(1,1,1,1))
  
  # tinnitus
  tinnitus <- subset(gwasResults, select = c("chr","rsid","pos","tinnitus.log10p"))
  gwasResults_allChr__tinnitus <- rbind(gwasResults_allChr__tinnitus,tinnitus)

  par(old.par)
  dev.off()
}

# GENOME-WIDE: tinnitus
jpeg(file="tinnitus_QQPLOT", width=2000,height=1000)
plotPvals("tinnitus",gwasResults_allChr__tinnitus,TRUE,FALSE)
dev.off()
jpeg(file="tinnitus_MANHATTAN", width=2000,height=1000)
plotPvals("tinnitus",gwasResults_allChr__tinnitus,FALSE,TRUE)
dev.off()
