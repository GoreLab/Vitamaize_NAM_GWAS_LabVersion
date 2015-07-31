rm(list = ls())

###Required files:
### (1) QTL array file
### (2) Imputed genotype matrices from get.me.my.snps.4 perl (GMMS4) script
### (3) GWAS summary files (Complete.Results...)

dir.for.all.imputed.SNPs <-"C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_GMMS4_imputed_matrix_T3_2015\\"
dir.for.GWAS.results <- "C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(10)GWAS Analysis\\RUV GWAS 25fam_alldata_alpha01_2015_FINAL\\"

#generate list of QTL numbers and common support intervals for GMMS4 script
all.JL.QTL.SI <- read.table("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Common_SI_array_for_tri_auto_June_2015.txt", head=TRUE)

for (k in c(9)){         #where k is total number of common support intervals  , large files for k=c(9,10,15,19,20,39)
  QTL <- all.JL.QTL.SI[k,1]
  gene <- as.character(all.JL.QTL.SI[k,5])
  chr <- all.JL.QTL.SI[k,2]
  left.bound <- all.JL.QTL.SI[k,3]
  right.bound <- all.JL.QTL.SI[k,4]
  common.SI <- paste(left.bound,"-",right.bound, sep = "")

  trait.list <- NULL
  for (m in 6:13){    #columns with significant traits
    select.trait <- as.character(all.JL.QTL.SI[k,m])
    trait.list <- c(trait.list, select.trait)
    #if(select.trait != "NA"){trait.list <- c(trait.list, select.trait)} else {break}     ### didn't work
    #while(select.trait != NA){trait.list <- c(trait.list, select.trait) ; if (select.trait == "NA") break}   ### didn't work
    
  } #end m traits
  trait.list <- trait.list[!is.na(trait.list)]
  
  #create a directory for each common QTL
dir.create(paste("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_GMMS4_imputed_matrix_T3_2015\\QTL_",QTL,"_imputed_matrix_for_",gene, sep=""))
output.folder <- (paste("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_GMMS4_imputed_matrix_T3_2015\\QTL_",QTL,"_imputed_matrix_for_",gene, "\\", sep=""))
for (i in trait.list){
     
   #identify row with corresponding QTL and trait under investigation to obtain individual QTL support interval
   row.of.interest <- tabular.summary[which((tabular.summary[,14] == QTL) & (tabular.summary[,1]) == i),]
   indiv.QTL.left.bound <- row.of.interest[1,4]
   indiv.QTL.right.bound <- row.of.interest[1,5]
   indiv.QTL.SI <- paste(indiv.QTL.left.bound,"-",indiv.QTL.right.bound, sep = "")
  
  #read in imputed genotype matrix file for specified chromosomal region from GMMS4
    #NOTE for GMMS4 this selections for individual QTL support interval
  setwd(dir.for.all.imputed.SNPs)
  GWAS.SNPs.imputed.all <- read.delim(paste("results_file_QTL",QTL,"_trait_",i,"_chr_",chr,"_region_",indiv.QTL.SI,".txt", sep=""))
  
  #####exception - file for chromosome 9
   ####GWAS.SNPs.imputed.all <- read.delim("results_file_QTL9_chr_2_region_54090079-80090079.txt")
  
  #read in sig RMIP GWAS SNPs by trait
  setwd(paste(dir.for.GWAS.results,"\\",i,"_iterations\\",sep = ""))
  GWAS.results <- read.delim(paste("Complete_Results_",i,"_with_effect_estimates.txt", sep=""))
    
    #obtain only sig RMIP SNPs
    #GWAS.results <- GWAS.results[-which(GWAS.results[,ncol(GWAS.results)] < 5),]      #applicable when last column of RMIP file is the RMIP value
    GWAS.results <- GWAS.results[-which(GWAS.results[,4] < 5),]
    
    #sort by chromosome and bp
    GWAS.results <- GWAS.results[order(GWAS.results[,2]),]
    GWAS.results <- GWAS.results[order(GWAS.results[,1]),]
  
    #obtain results only for chromosome of interest
    GWAS.results <- GWAS.results[which(GWAS.results[,1] == chr),]
    query.SNPs <- as.matrix(GWAS.results[,2])
   
   #reduce number of SNPs to be investigated from those in common support interval (in GWAS.SNPs.imputed.RMIPsig) to those in individual support interval (corresponding to number of SNPs in tab summary) 
    SNP.list.indiv.SI <- NULL
    for (m in 1:nrow(query.SNPs)){
      if((query.SNPs[m,1] > indiv.QTL.left.bound) & (query.SNPs[m,1] < indiv.QTL.right.bound)){
      SNP.list.indiv.SI <- rbind(SNP.list.indiv.SI, query.SNPs[m,])}  
     } # end m rows  
  
  if(nrow(SNP.list.indiv.SI) < 1){break
    }else{
  
  GWAS.SNPs.imputed.all <- as.matrix(GWAS.SNPs.imputed.all)
  SNP.list.indiv.SI <- as.matrix(SNP.list.indiv.SI)

  SNP.list <- NULL
  GWAS.SNPs.imputed.RMIPsig <- NULL
    for (j in 1:nrow(SNP.list.indiv.SI)){
      #obtain imputed genotype matrices for significant SNPs
       SNP.list <- GWAS.SNPs.imputed.all[which(GWAS.SNPs.imputed.all[,4]==SNP.list.indiv.SI[j,1]),]
       GWAS.SNPs.imputed.RMIPsig <-  rbind(GWAS.SNPs.imputed.RMIPsig, SNP.list)
      } #end j in nonimputed rows

    #output SNPs within individual support interval; note that this number may still be higher than that in tab summary if there are multiple allelic states for a given SNP
    write.table(GWAS.SNPs.imputed.RMIPsig, paste(output.folder,"imputed_RMIP_SNPs_",gene,"_interval_",i,"_QTL",QTL,"_chr",chr,".txt", sep=""), sep = "\t", row.names = FALSE)
 
 } #end if/else loop accounting for zero SNPs in region
} #end i in trait.list

} #end k in QTL sequence

