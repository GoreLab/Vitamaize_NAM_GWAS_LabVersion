rm(list = ls())

###Required files:
### (1) QTL array file
### (2) Imputed genotype matrices from get.me.my.snps.4 perl (GMMS4) script
### (3) GWAS summary files (Complete.Results...)

dir.for.all.imputed.SNPs <-"C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis\\Results_from_GMMS4_imputed_matrix_T3_2015\\Files_from_GMMS4\\"
dir.for.GWAS.results <- "C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(10)GWAS Analysis\\RUV GWAS 25fam_HMPonly_TASSEL3_alpha01_2015 corr\\"

#generate list of QTL numbers and common support intervals for GMMS4 script
all.JL.QTL.SI <- read.table("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis\\Common_SI_array_for_tri_auto_T3_2015_SI01_new_GWAS.txt", head=TRUE)

for (k in 1:38){         #where k is total number of common support intervals
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
  

#read in imputed genotype matrix file for specified chromosomal region from GMMS4
    #NOTE for GMMS4 select common JL support interval for region of interest

setwd(dir.for.all.imputed.SNPs)
GWAS.SNPs.imputed.all <- read.delim(paste("results_file_QTL",QTL,"_chr_",chr,"_region_",common.SI,".txt", sep=""))

#read in sig RMIP GWAS SNPs by trait

dir.create(paste("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis\\Results_from_GMMS4_imputed_matrix_T3_2015\\QTL_",QTL,"_imputed_matrix_for_",gene, sep=""))
output.folder <- (paste("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis\\Results_from_GMMS4_imputed_matrix_T3_2015\\QTL_",QTL,"_imputed_matrix_for_",gene, "\\", sep=""))
for (i in trait.list){
  setwd(paste(dir.for.GWAS.results,"\\",i,"\\",sep = ""))
                                    
    GWAS.results <- read.delim(paste("Complete_Results_",i,".txt", sep=""))
    
    #obtain only sig RMIP SNPs
    GWAS.results <- GWAS.results[-which(GWAS.results[,ncol(GWAS.results)] < 5),]
    
    #sort by chromosome and bp
    GWAS.results <- GWAS.results[order(GWAS.results[,2]),]
    GWAS.results <- GWAS.results[order(GWAS.results[,1]),]
  
    #obtain results only for chromosome of interest
    GWAS.results <- GWAS.results[which(GWAS.results[,1] == chr),]
  
  
    ###NOTE!!! In some cases there will not be any non-imputed RMIP SNPs in the following list, causing the script to stall
    ###A try function may need to be instituted to make this flow seamlessly
    ###Currently, upon obtaining error, Trait QTL is checked for # of GWAS hits in interval from master file; if zero then restart loop with next QTL
    ###Such was the case with carot QTL9 and 10
  #GWAS.SNPs.nonimputed.RMIPsig <- read.delim(paste("Output_GWAS_Results_",i,"_lt_0.05_Chr_",chr,".txt",sep=""), head=TRUE)
  query.SNPs <- as.matrix(GWAS.results[,2])
  GWAS.SNPs.imputed.all <- as.matrix(GWAS.SNPs.imputed.all)

  SNP.list <- NULL
  GWAS.SNPs.imputed.RMIPsig <- NULL
    for (j in 1:nrow(query.SNPs)){
      #obtain SNPs from column 1 in GWAS.SNPs.nonimputed.RMIPsig that are in GWAS.SNPs.imputed.all
       SNP.list <- GWAS.SNPs.imputed.all[which(GWAS.SNPs.imputed.all[,4]==query.SNPs[j,1]),]
       #SNP.list <- GWAS.SNPs.imputed.all[which(GWAS.SNPs.imputed.all[,1] % GWAS.SNPs.nonimputed.RMIPsig[,1])]
       #SNP.list <- subset(GWAS.SNPs.imputed.all, GWAS.SNPs.imputed.all[,1] == GWAS.SNPs.nonimputed.RMIPsig[,1])
      GWAS.SNPs.imputed.RMIPsig <-  rbind(GWAS.SNPs.imputed.RMIPsig, SNP.list)
      } #end j in nonimputed rows

    GWAS.SNPs.imputed.RMIPsig.unique  <- unique(GWAS.SNPs.imputed.RMIPsig)     #this pulls the exact SNPs from which the GWAS list was derived - without unique function, duplicates are generated
    write.table(GWAS.SNPs.imputed.RMIPsig.unique, paste(output.folder,"imputed_RMIP_SNPs_",gene,"_interval_",i,"_QTL",QTL,"_chr",chr,"REDO.txt", sep=""), sep = "\t", row.names = FALSE)

    ############NOTE: script currently needs modification to allow tables to be exported for entire vector of traits
    ############Currently, you can only run it one trait.list at a time because of some info that needs to be cleared from memory
        ###############before reentering the loop

} #end i in trait.list

} #end k in QTL sequence

