rm(list = ls())

###########NOTE this has been changed from original on 04/01/15 to only contain correlation analysis for JL markers
### Follow this analysis up by selecting which markers appear to be MC and excluding them from marker vector which will be used in TASSEL rescan function
###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version) and residual row removed from file

setwd("C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/JL/") #changed CHD 5-11
home.dir <- getwd()
geno.path = paste(home.dir, "/validate_CBK.AEL/",sep='')
#trait <- "TRAIT5_ACAR_RUV" #changed CHD 5-11--commented out
trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs") #added CHD 5-11
BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"
testing.correlation <- TRUE
transformed = FALSE

if(transformed == FALSE){
  #For JL on untransformed BLUEs, only for backtransformation validation
  JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/JL_Analysis_Scripts.Results/UNTRANSF_JLResults/",sep='')
  correlation.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/JL_Analysis_Scripts.Results/UNTRANSF_JLResults/",sep='')
  pheno.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Validation_Backtransformation/BLUES_Untransf/",sep='')
}else{
#For master JL analysis
  JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Final_Models_forMultiCollCheck/",sep='')
  correlation.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Correlation_Matrices/",sep='')
  pheno.path = paste(home.dir, "/Phenotypes/",sep ='')
}

library(multtest)

#Read in the gneotypic data
#setwd(paste(home.dir, "\\validate_CBK.AEL", sep = "")) #changed CHD 5-11
genotypes <- read.table(paste(geno.path,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Eventually, loop around the traits. This loop will begin here

for(i in trait){

#Read in the trait under study
  if(transformed == FALSE){
    pheno.data = read.table(paste(pheno.path,"NAM_Carot_BestTraits_BLUEs_Merged_All_Env_",i,"_RUV.txt",sep=''),head=TRUE)
  }else{
    pheno.data <- read.table(paste(pheno.path,"BLUEs_No_Outliers_Transformed_all_for_TASSEL_",i,".txt", sep = ''), head = TRUE) #changed CHD file name 5-11 to match toco traits
  }
head(pheno.data)

#Read in the results from JL analysis
if(transformed == FALSE){
  TASSEL.model.results = read.table(paste(JL.path,i,"_model_3fromSF_0.01_UNTRANSF_Final_forR.txt",sep=''),head=TRUE)
}else{
#setwd(paste(home.dir, "\\validate_CBK.AEL\\CHD_Tassel3fromSF_modified0.01\\R.JL.no.MultiColl\\Toco_Models_postMCCorr-whenApplicable", sep = "")) #changed CHD to direct path from top of script 5-12
#new path for SI01 required here 
  TASSEL.model.results <- read.table(paste(JL.path,i,"_model_3fromSF_0.01_Final_R.formatted.txt", sep = "") , head = TRUE)
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
#Fit the JL model and obtain pop*marker terms
#attach(geno.and.pheno) #CHD commented out on 5-12; not needed and was resulting in masking
#setwd(paste(home.dir, "\\JL_Results_", BLUE.or.BLUP,"\\Pop_by_Marker_Terms_from_R\\",  sep = ""))
 

#Calculate a correlation matrix between correlated SNPs, and print them out to a data set.
if(testing.correlation == TRUE){
  #setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(19) Evaluating Multicollinearity\\")
  #Description.of.SNPs <- "Chr_ALL_SNPs_ACAR"
  #SNP.IDs <- c("S_44638563", "S_48526286", "S_118969928", "S_136572666")
  #uncomment two lines below if testing correlation across all significant markers
  TASSEL.model.results.2 <- TASSEL.model.results[2:nrow(TASSEL.model.results),]
  SNP.IDs <- as.vector(unique(TASSEL.model.results.2[,4]))
  Description.of.SNPs = "allSignifSNPs" #line added by CHD 5-12
  #Description.of.SNPs <- "Chr_5_SNPs_bcar"
  #SNP.IDs <- c("S_57225641", "S_56737096", "S_59271817")
  #Description.of.SNPs <- "Chr_1_SNPs_totcar"
  #SNP.IDs <- c("S_289019777", "S_287963456")
  
  SNP.set <- geno.and.pheno[,which(colnames(geno.and.pheno) %in% SNP.IDs)]
  correlation.matrix <- cor(SNP.set, use = "pairwise.complete.obs")
  write.table(correlation.matrix, paste(correlation.path,"Correlation.matrix_",Description.of.SNPs,"_",i,".txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
}

#detach(geno.and.pheno) #CHD added then immediately commented out on 5-12; should be used if using attach(), but better to avoid entirely.

}# End for(i in trait)