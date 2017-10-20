setwd("/Users/anybody/Documents/Toco_NAM/")
home.dir <- getwd()
GWAS.results.dir = paste(home.dir,"/GWAS/Summaries_by_Trait/",sep='')
tabSummary.path = "/Summary_Tables.Figures/"

all.signif.RMIP.SNPs = read.table(paste(GWAS.results.dir,"Complete_Results_All.Traits_RMIPgt4_dT3.Redone.txt",sep=''),head=TRUE,stringsAsFactors=FALSE)
setwd(paste(home.dir,tabSummary.path,sep=''))
commonSI.summ = read.table("Tab_Sum_Final_dT3_Redone_with_Common_SI_Info_left_bound.txt",stringsAsFactors=FALSE)

which.cSI.this.RMIP.vect = rep(NA,nrow(all.signif.RMIP.SNPs))
for(SNP in (1:nrow(all.signif.RMIP.SNPs))){
  this.row = all.signif.RMIP.SNPs[SNP,]
  chr.SNP = as.character(this.row[2])
  trait.SNP = as.character(this.row[1])
  pos.SNP = as.numeric(this.row[3])
  
  RMIP.SNP.vect = c(trait.SNP,chr.SNP)
  #if chr,pos,trait of SNP matches that combo for any row in master.commonSI.summary 
  trait.chr.match = subset(commonSI.summ,commonSI.summ[,1] == trait.SNP & commonSI.summ[,2] == chr.SNP)
  print(nrow(trait.chr.match))
  if(nrow(trait.chr.match)==0){
    next
  } else {
  for (match.row in (1:nrow(trait.chr.match))){
    this.row.cSI = trait.chr.match[match.row,]
      indiv.trait.SI.left = as.numeric(this.row.cSI[4])
      indiv.trait.SI.right = as.numeric(this.row.cSI[5])
    if (pos.SNP >= indiv.trait.SI.left && pos.SNP <= indiv.trait.SI.right){
      which.cSI.this.RMIP.vect[SNP] = as.numeric(this.row.cSI[12])
      print(pos.SNP)
      print(this.row.cSI[12])
    }
  }
}
}

RMIP.annot.by.cSI = cbind(all.signif.RMIP.SNPs,which.cSI.this.RMIP.vect)
write.table(RMIP.annot.by.cSI,paste(GWAS.results.dir,"Complete_Results_allTraits_RMIPgt4_dT3.Redone_with_QTLnumber.txt",sep=''),sep='\t',row.names=FALSE,quote=FALSE)