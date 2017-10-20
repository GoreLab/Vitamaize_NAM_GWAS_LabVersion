rm(list = ls())

###Required files:
### (1) QTL array file
### (2) Imputed genotype matrices from get.me.my.snps.3.1 perl (GMMS3.1) script
### (3) GWAS summary files (Complete.Results...)
### (4) Tabular summary file containing individual QTL support intervals

setwd("/Users/anybody/Documents/Toco_NAM/")
home.dir <- getwd()
expression.dir = paste(home.dir,"/Expression/",sep='')
dir.for.5.0.SNPs <-paste(home.dir,"/GWAS/Toco_GWAS_SNPs_dT3.Redone/",sep='')
dir.for.5.1.matrices = paste(expression.dir,"ImputedMatricesfrom5.1/",sep='')
dir.for.GWAS.results <- paste(home.dir,"/GWAS/Summaries_by_Trait/",sep='')
tabSummary.path = paste(home.dir,"/Summary_Tables.Figures/",sep='')
numb.common.SI = 52
max.numb.traits.per.SI = 10

#generate list of QTL numbers and common support intervals for GMMS4 script
all.JL.QTL.SI <- read.table(paste(tabSummary.path,"Common_SI_array_for_tri_auto_June_2016.txt",sep=''), head=TRUE,stringsAsFactors=FALSE)

SI.numbs = seq(1,numb.common.SI)
#SI.numbs.needing.Redo = c(1,2,4,5,8,12,14,18,21,23,24,26,30,33,36,37,38,41,43,45,48,49,52) #used for dT3 redo

for (k in SI.numbs){ #changed this to SI.numbs.needing.Redo for dT3 redo
  QTL <- all.JL.QTL.SI[k,1]
  gene <- as.character(all.JL.QTL.SI[k,5])
  chr <- all.JL.QTL.SI[k,2]
  left.bound <- all.JL.QTL.SI[k,3]
  right.bound <- all.JL.QTL.SI[k,4]
  common.SI <- paste(left.bound,"-",right.bound, sep = "")

  trait.list <- NULL
  for (m in 6:(5+max.numb.traits.per.SI)){      #added by CHD since one toco interval had 9 associated traits; first trait column is column 6 (index 6)
    select.trait <- as.character(all.JL.QTL.SI[k,m])
    trait.list <- c(trait.list, select.trait)   
  } #end m traits
  trait.list <- trait.list[!is.na(trait.list)]
  
#read in imputed genotype matrix file for specified chromosomal region from GMMS4
    #NOTE for GMMS4 select common JL support interval for region of interest

setwd(paste(dir.for.5.0.SNPs,sep=''))
GWAS.SNPs.imputed.all <- read.delim(paste("results_file_QTL",QTL,"_chr_",chr,"_region_",common.SI,".txt", sep=""),sep="\t",head=T,stringsAsFactors=FALSE,strip.white=TRUE)
print(as.character(GWAS.SNPs.imputed.all[1,4]))
print(dim(GWAS.SNPs.imputed.all))
#read in sig RMIP GWAS SNPs by trait

dir.create(paste(dir.for.5.1.matrices,"QTL_",QTL,"_imputed_matrix_for_chr",chr, sep=""))
output.folder <- paste(dir.for.5.1.matrices,"QTL_",QTL,"_imputed_matrix_for_chr",chr, sep="")
for (i in trait.list){
  
  #to identify SNPs specifically within common support interval read in tab summary
  tabular.summary <- read.table(paste(tabSummary.path,"Tab_Sum_Final_dT3_Redone_with_Common_SI_Info_left_bound.txt",sep=''), head = TRUE,stringsAsFactors=FALSE)
  
  #identify row with corresponding QTL and trait under investigation to obtain individual QTL support interval
  row.of.interest <- tabular.summary[which((tabular.summary[,12] == QTL) & (tabular.summary[,1]) == i),]
  indiv.QTL.left.bound <- row.of.interest[1,4]
  indiv.QTL.right.bound <- row.of.interest[1,5]
  
  if(i == "Total_Tocopherols"){i = "totalT"}
  if(i == "Total_Tocotrienols"){i = "totalT3"}
  if(i == "Total_Tocochromanols"){i = "totalTocochrs"}
  setwd(paste(dir.for.GWAS.results,sep = ""))
  
  if(i == "dT3"){
      GWAS.results <- read.delim(paste("Complete_Results_",i,"_Redone.txt", sep=""),sep="\t")
  }else{
    GWAS.results <- read.delim(paste("Complete_Results_",i,".txt", sep=""),sep="\t")
  }
  
    #obtain only sig RMIP SNPs
    #GWAS.results <- GWAS.results[-which(GWAS.results[,ncol(GWAS.results)] < 5),]      #applicable when last column of RMIP file is the RMIP value
    GWAS.results <- GWAS.results[-which(GWAS.results[,4] < 5),]
    
    #sort by chromosome and bp
    GWAS.results <- GWAS.results[order(GWAS.results[,2]),]
    GWAS.results <- GWAS.results[order(GWAS.results[,1]),]
  
    #obtain results only for chromosome of interest
    GWAS.results <- GWAS.results[which(GWAS.results[,1] == chr),]
    query.SNPs <- as.matrix(GWAS.results[,2])
  
    ###NOTE!!! In some cases there will not be any non-imputed RMIP SNPs in the following list, causing the script to stall
    ###A try function may need to be instituted to make this flow seamlessly
    ###Currently, upon obtaining error, Trait QTL is checked for # of GWAS hits in interval from master file; if zero then restart loop with next QTL
  
   #reduce number of SNPs to be investigated from those in common support interval (in GWAS.SNPs.imputed.RMIPsig) to those in individual support interval (corresponding to number of SNPs in tab summary) 
    SNP.list.indiv.SI <- NULL
    for (m in 1:nrow(query.SNPs)){
      if((query.SNPs[m,1] > indiv.QTL.left.bound) & (query.SNPs[m,1] < indiv.QTL.right.bound)){
      SNP.list.indiv.SI <- rbind(SNP.list.indiv.SI, query.SNPs[m,])}  
     } # end m rows  
  
  if(is.null(nrow(SNP.list.indiv.SI))){
    print(paste("No snps in interval ",common.SI," for trait ",i,sep=''))
    next
  }else{
  
  GWAS.SNPs.imputed.all <- as.matrix(GWAS.SNPs.imputed.all)
  SNP.list.indiv.SI <- as.matrix(SNP.list.indiv.SI)

  SNP.list <- NULL
  GWAS.SNPs.imputed.RMIPsig <- NULL
    for (j in 1:nrow(SNP.list.indiv.SI)){
      print(SNP.list.indiv.SI[j,1])
      #obtain imputed genotype matrices for significant SNPs
       SNP.list <- GWAS.SNPs.imputed.all[which(GWAS.SNPs.imputed.all[,4]==as.character(SNP.list.indiv.SI[j,1])),]
       print(dim(SNP.list))
       GWAS.SNPs.imputed.RMIPsig <-  rbind(GWAS.SNPs.imputed.RMIPsig, SNP.list)
      } #end j in nonimputed rows

    #output SNPs within individual support interval; note that this number may still be higher than that in tab summary if there are multiple allelic states for a given SNP
    write.table(GWAS.SNPs.imputed.RMIPsig, paste(output.folder,"/imputed_RMIP_SNPs_interval_",i,"_QTL",QTL,"_chr",chr,".txt", sep=""), sep = "\t", row.names = FALSE)
  
  } #end if/else loop accounting for zero SNPs in region
} #end i in trait.list

} #end k in QTL sequence