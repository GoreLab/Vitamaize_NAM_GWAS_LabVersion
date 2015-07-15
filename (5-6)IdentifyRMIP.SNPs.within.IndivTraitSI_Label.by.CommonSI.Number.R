

setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
home.dir <- getwd()
GWAS.results.dir <- "/GWAS_Analysis/GWAS_25fam_HMPonly_TASSEL3_alpha01_2015_corr/"
tabSummary.path = "/Tabular_Summaries/"

setwd(paste(home.dir,GWAS.results.dir,sep=''))
all.signif.RMIP.SNPs = read.table("Complete_Results_allTraits_RMIPgt4_forR.txt",as.is=TRUE,head=TRUE)

setwd(paste(home.dir,tabSummary.path,sep=''))
commonSI.summ = read.table("Tab_Sum_tocos_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions_LODscores.txt",stringsAsFactors=FALSE)

which.cSI.this.RMIP.vect = rep(NA,nrow(all.signif.RMIP.SNPs))
for(SNP in (1:nrow(all.signif.RMIP.SNPs))){
  this.row = all.signif.RMIP.SNPs[SNP,]
  chr.SNP = as.character(this.row[1])
  trait.SNP = as.character(this.row[5])
  if(trait.SNP == "totalT"){trait.SNP = "Total_Tocopherols"}
  pos.SNP = as.numeric(this.row[2])
  
  RMIP.SNP.vect = c(trait.SNP,chr.SNP)
  #if chr,pos,trait of SNP matches that combo for any row in master.commonSI.summary 
  trait.chr.match = subset(commonSI.summ,commonSI.summ[,1] == trait.SNP & commonSI.summ[,2] == chr.SNP)
  print(nrow(trait.chr.match))
  if(nrow(trait.chr.match)==0){
    next
  } else {
  for (match.row in (1:nrow(trait.chr.match))){
    this.row.cSI = trait.chr.match[match.row,]
      cSI.left = as.numeric(this.row.cSI[4])
      cSI.right = as.numeric(this.row.cSI[5])
    if (pos.SNP >= cSI.left && pos.SNP <= cSI.right){
      which.cSI.this.RMIP.vect[SNP] = as.numeric(this.row.cSI[12])
      print(pos.SNP)
      print(this.row.cSI[12])
    }
  }
}
}

RMIP.annot.by.cSI = cbind(all.signif.RMIP.SNPs,which.cSI.this.RMIP.vect)
write.table(RMIP.annot.by.cSI,paste(home.dir,tabSummary.path,"Complete_Results_allTraits_RMIPgt4_with_QTLnumber.txt",sep=''),sep='\t')