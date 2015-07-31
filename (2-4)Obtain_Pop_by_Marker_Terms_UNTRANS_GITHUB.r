rm(list = ls())

###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version) and residual row removed from file



#Use the multtest library to obtain FDR-adjusted P-values


setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")                                #changed CHD 5-11
home.dir <- getwd()

trait.set = "Carot" #options: "Carot" or "Toco"

geno.path = "C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/JL/validate_CBK.AEL/"
UNTRANS.pheno.path = paste(home.dir, "/Validation_Backtransformation/BLUES_Untransf/",sep ='')
JL.path = paste(home.dir,"/R.JL.no.MultiColl/",trait.set,"_FinalModels_postMCCorrPostRescan-whenApplicable/",sep='')
popByMarker.path = paste(home.dir, "/Validation_Backtransformation/Untrans_Allelic.Effect.Estimates/", sep = "")

library(multtest)

#Read in the genotypic data
#setwd(paste(home.dir, "\\validate_CBK.AEL", sep = ""))                                       #CHD commented out 5-14
genotypes <- read.table(paste(geno.path,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Eventually, loop around the traits. This loop will begin here                                #CHD 5-14 loop made so that carots or tocos can be run automatically from top of script.
if (trait.set == "Carot"){
  trait <- c("LUT","ZEI","BCRY","ACAR","PHYF","THLYC","ZEA", "BCAR", "TOTCAR")    
  trait.collinear <- c("ZEA", "BCAR", "TOTCAR") 
} else if (trait.set == "Toco"){
  trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs")       #CHD added 5-11
  trait.collinear <- c("aT3","dT3","gT3","totalT3","totalTocochrs")   
}else{print("Error: no trait set selected. Please specify carots or tocos.")}
  
#BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"

for(i in trait){
#Read in the trait under study
#setwd(paste(home.dir, "\\Phenotypes\\",sep = ""))
  if (trait.set == "Carot"){
      pheno.data <- read.table(paste(UNTRANS.pheno.path,"NAM_Carot_BestTraits_BLUEs_Merged_All_Env_",i,"_RUV.txt", sep = ""), head = TRUE)
  } else if (trait.set == "Toco"){
      pheno.data <- read.table(paste(UNTRANS.pheno.path,"BLUEs_Untransformed_No_Out_201304019_for_R_",i,".txt", sep = ""), head = TRUE)
  }

head(pheno.data)

#Read in the results from JL analysis
#setwd(paste(home.dir, "\\Data_for_1pct_Support_Intervals\\JL_output\\JL_model_modified_SI01\\", sep = ""))
#new path for SI01 required here (CBK note)
if (trait.set == "Carot"){
  if (i %in% trait.collinear){
    TASSEL.model.results <- read.table(paste(JL.path,"MC_corrected_",i,"_forR.txt", sep = "") , head = TRUE)
  } else {
    TASSEL.model.results <- read.table(paste(JL.path,i,"_finalModel_withModelR2_forR.txt", sep = "") , head = TRUE)
  }
}
  
if (trait.set == "Toco"){
  if (i %in% trait.collinear){   
    TASSEL.model.results <- read.table(paste(JL.path,"MC_corrected_",i,"_postRescan_R.formatted.txt", sep = "") , head = TRUE)        #CHD added loop to handle both MC and non-MC traits
  } else {
    TASSEL.model.results <- read.table(paste(JL.path,i,"_model_3fromSF_0.01_Final_R.formatted.txt", sep = "") , head = TRUE)
  }
}

#### #########################################################################################################################################
#Get the appropriate SNPs and merge the phenotypic and genotypic data together.
  
#Parse out all of the SNPs that were selected in the model
geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)]
geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))



#pheno.data will always have more data becuase IBM is included in the phenotypic data.
colnames(pheno.data)[1] = "Geno_Code" #line added by CHD 5-14; toco pheno files are formatted for TASSEL so use <trait> header format; this line will unify with carot format for the column needed for merge.
geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")
print(dim(geno.and.pheno))

#Add a population column
geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
colnames(geno.and.pheno)[2] <- "pop"


##############################################################################################################################################
#Fit the JL model and obtain pop*marker terms
#attach(geno.and.pheno)   #cHD commented out 5/14; attaching causes masking and other problems, avoiding. geno.and.pheno was already specified below where needed.


#Get the model order of the model in R to be the same as that in TASSEL.
seq <- (seq(1:nrow(TASSEL.model.results))-1)[-1]                 #BY CBK: obtain number of sig SNPs from JL model, but remove 1st pop row) 
model.order <- cbind(seq, TASSEL.model.results[-1,c(2,3,4)])     #BY CBK: bind num of sig SNPs with model results (chr, position, SNP_name) 
#Sort by chromosome and bp
model.order <- model.order[order(model.order[,3]),]              #BY CBK: order model based on position (seq now different)
model.order <- model.order[order(model.order[,2]),]              #BY CBK: order model by chromosome (and position)
model.order <- cbind(model.order, seq(1:nrow(model.order)))      #BY CBK: add new sequence numbering to existing order model 
#Sort by the first column of marker order so that the markers will be put into the model in the same order as they were selected.
model.order <- model.order[order(model.order[,1]),]              #BY CBK: obtain original numbering
index <- model.order[,5]+3 #3 is added so because the first SNP is on the fourth column             #by CBK: just the model order by chr and position

if(i=="totalT"){i="Total_Tocopherols"}
if(i=="totalT3"){i="Total_Tocotrienols"}
if(i=="totalTocochrs"){i="Total_Tocochromanols"}
base.model <- paste(i, "_RUV ~ pop",sep = "")                    #By CBK: paste trait and pop+   - note that tilde means "modeled as"; CHD removed "+" after pop 5-14 because + already appended with each marker
for(k in index){                                                 #BY CBK: cycle through rows in index, which is the sequence ordered by chr and position
  base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")   #BY CBK:  writing out model of trait modeled by pop plus pop by SNPs in model
} #end for(k in 4:ncol(geno.and.pheno))

#JL.model <- lm(AT ~ pop+pop:S_38836669, data = geno.and.pheno) 
JL.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno)    # linear model of trait modeled by pop, pop:SNPs specifying data in geno.and.pheno
print(i)
print(summary(JL.model)$df)
print(summary(JL.model)$na.action)

pop.by.marker.effects <- summary(JL.model)$coefficients[-c(1:25),]      # obtain descriptive statistics for linear model and model coefficients for pop and SNP within pop, access individual columns with model coefficients
         #CBK changed 26 back to 25 - THIS CONTRIBUTED TO THE LOSS OF POP1:snp1 EFFECT IN ALL POP.BY.MKR.EFFECT FILES
head(pop.by.marker.effects)                                              #header for this data matrix pop term, estimate, stderr, tvalue, prob, FDR adj pvalue

#Perform the FDR procedure on the P-values from pop.by.marker.effects
 res <- mt.rawp2adjp(pop.by.marker.effects[,4], "BH")
 adjp <- res$adjp[order(res$index), ]
                                       
 pop.by.marker.effects <- cbind(pop.by.marker.effects, adjp[,2])
 colnames(pop.by.marker.effects)[5] <- "FDR_Adjusted_P-values"

#Write out the results, put the in the same directory as the JL results from TASSEL
 pop.by.marker.effects <- cbind(rownames(pop.by.marker.effects), pop.by.marker.effects)
 colnames(pop.by.marker.effects)[1] <- "Term"
 
if(i=="Total_Tocopherols"){i="totalT"}
if(i=="Total_Tocotrienols"){i="totalT3"}
if(i=="Total_Tocochromanols"){i="totalTocochrs"}
 write.table(pop.by.marker.effects, paste(popByMarker.path,"UNTRANS_Pop.by.Marker.Effect.Estimates.from.R.",i,".SI01.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

}# End for(i in trait)