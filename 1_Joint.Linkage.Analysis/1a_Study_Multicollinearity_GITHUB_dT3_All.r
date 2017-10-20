rm(list = ls())

###########NOTE this has been changed from original on 04/01/15 to only contain correlation analysis for JL markers
### Follow this analysis up by selecting which markers appear to be MC and excluding them from marker vector which will be used in TASSEL rescan function
###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version) and residual row removed from file

setwd("/Users/anybody/Documents/Toco_NAM/") #changed CHD 5-11
home.dir <- getwd()
geno.path = paste(home.dir, "/Geno.Pheno_Inputs/",sep='')
trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs") #added CHD 5-11
#trait = "dT3" #CHD added 19 May for re-do
dT3.new.path = paste(home.dir,"/Methods/dT3_removeExtremeVal_test/new.trans_new.perm_FINAL/",sep='')
correlation.path = dT3.new.path
BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"
testing.correlation <- TRUE
transformed = TRUE

if(transformed == FALSE){
  #For JL on untransformed BLUEs, only for backtransformation validation
  JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/JL_Analysis_Scripts.Results/UNTRANSF_JLResults/",sep='')
  pheno.path = "/Users/anybody/Documents/Toco_NAM/Geno.Pheno_Inputs/BLUES_Untransf/"
  #pheno.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Validation_Backtransformation/BLUES_Untransf/",sep='')
}else{
#For master JL analysis
  #JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Final_Models_forMultiCollCheck/",sep='')
  JL.path = paste(home.dir,"/JL/Toco_FinalModels_postMCCorrPostRescan/",sep='')
  #pheno.path = paste(home.dir, "/Phenotypes/",sep ='')
  pheno.path = paste(home.dir, "/Geno.Pheno_Inputs/",sep='')
}

library(multtest)

#Read in the gneotypic data
#setwd(paste(home.dir, "\\validate_CBK.AEL", sep = "")) #changed CHD 5-11
genotypes <- read.table(paste(geno.path,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Eventually, loop around the traits. This loop will begin here

for(i in trait){

#Read in the trait under study
  if(i == "dT3"){
    pheno.data = read.table(paste(dT3.pheno.path,"BLUEs_No_Outliers_Transformed_dT3.only.Redone.txt",sep=''),head=TRUE,stringsAsFactors=FALSE)
    pheno.data$dT3 = as.numeric(pheno.data$dT3)
  }else{
    pheno.data = read.table(paste(pheno.path,"BLUEs_No_Outliers_Transformed_all_for_TASSEL_",i,".txt", sep = ""), head = TRUE)
  }
  
  head(pheno.data)

#Read in the results from JL analysis
if(transformed == FALSE){
  TASSEL.model.results = read.table(paste(JL.path,i,"_model_3fromSF_0.01_UNTRANSF_Final_forR.txt",sep=''),head=TRUE)
}else{
  if(i == "dT3"){
    TASSEL.model.results <- read.table(paste(JL.path,"JL_fromTrans_dT3_Redone_MC.corrected_Tassel3fromSF_forR.txt", sep = "") , head = TRUE,stringsAsFactors = FALSE)       
  }else if (i %in% trait.collinear){   
    print("dT3 went here...oops.")
    TASSEL.model.results <- read.table(paste(JL.path,"MC_corrected_",i,"_postRescan_R.formatted.txt", sep = "") , head = TRUE,stringsAsFactors = FALSE)        #CHD added loop 5-14 to handle both MC and non-MC traits
  } else {
    TASSEL.model.results <- read.table(paste(JL.path,i,"_model_3fromSF_0.01_Final_R.formatted.txt", sep = "") , head = TRUE,stringsAsFactors = FALSE)
  }
}

if(is.numeric(TASSEL.model.results[,2]) == FALSE){
 print(paste("The chromosome column in ", i," model results needs to be numeric. Please change this before proceeding.", sep = ""))
 break;
}

#### #########################################################################################################################################
#Get the appropriate SNPs and merge the phenotypic and genotypic data together.
  
#Parse out all of the SNPs that were selected in the model
geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)]
geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))

#pheno.data will always have more data becuase IBM is included in the phenotypic data.
colnames(pheno.data)[1] = "Geno_Code" #added by CHD 5-12 so that pheno data has header script is looking for (otherwise is "X.trait.")
geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

#Add a population column
geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
colnames(geno.and.pheno)[2] <- "pop"

##############################################################################################################################################
#Calculate a correlation matrix between correlated SNPs, and print them out to a data set.
if(testing.correlation == TRUE){
  #uncomment two lines below if testing correlation across all significant markers
  TASSEL.model.results.2 <- TASSEL.model.results[2:nrow(TASSEL.model.results),]
  SNP.IDs <- as.vector(unique(TASSEL.model.results.2[,4]))
  Description.of.SNPs = "allSignifSNPs" #line added by CHD 5-12
  
  SNP.set <- geno.and.pheno[,which(colnames(geno.and.pheno) %in% SNP.IDs)]
  correlation.matrix <- cor(SNP.set, use = "pairwise.complete.obs")
  write.table(correlation.matrix, paste(correlation.path,"Correlation.matrix_",Description.of.SNPs,"_",i,".txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
}
}# End for(i in trait)