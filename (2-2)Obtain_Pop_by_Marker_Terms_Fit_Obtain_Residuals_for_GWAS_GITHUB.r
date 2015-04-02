rm(list = ls())

###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version)



#Use the multtest library to obtain FDR-adjusted P-values


setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(9)JL Analysis\\Permutations")
home.dir <- getwd()


library(multtest)

#Read in the genotypic data
setwd(paste(home.dir, "\\GBS_SNPs", sep = ""))
genotypes <- read.table("imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt", head = TRUE)

#Eventually, loop around the traits. This loop will begin here

trait <- c("LUT_RUV", "ZEA_RUV","ZEI_RUV","BCRY_RUV","ACAR_RUV","PHYF_RUV","THLYC_RUV")    
trait.collinear <- c("BCAR_RUV", "TOTCAR_RUV")    

#BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"

for(i in trait){
#Read in the trait under study
setwd(paste(home.dir, "\\Phenotypes\\Phen_files_no_trait_prefix\\",sep = ""))
pheno.data <- read.table(paste(trait,"_BLUE_outtrans.txt", sep = ""), head = TRUE)
head(pheno.data)

#Read in the results from JL analysis
setwd(paste(home.dir, "\\Data_for_1pct_Support_Intervals\\JL_output\\JL_model_modified_SI01\\", sep = ""))
TASSEL.model.results <- read.table(paste(trait, "_JL_model_modified_SI01.txt"   , sep = "") , head = TRUE)

### if trait.collinear
#TASSEL.model.results <- read.table(paste(trait, "_JL_model_modified_no_multicollinearity_SI01.txt"   , sep = "") , head = TRUE)

#### #########################################################################################################################################
#Get the appropriate SNPs and merge the phenotypic and genotypic data together.
  
#Parse out all of the SNPs that were selected in the model
geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)]
geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))



#pheno.data will always have more data becuase IBM is included in the phenotypic data.

geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

#Add a population column
geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
colnames(geno.and.pheno)[2] <- "pop"


##############################################################################################################################################
#Fit the JL model and obtain pop*marker terms
attach(geno.and.pheno)


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

base.model <- paste(trait, " ~ pop+",sep = "")                    #By CBK: paste trait and pop+   - note that tilde means "modeled as"
for(k in index){                                                 #BY CBK: cycle through rows in index, which is the sequence ordered by chr and position
  base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")   #BY CBK:  writing out model of trait modeled by pop plus pop by SNPs in model
} #end for(k in 4:ncol(geno.and.pheno))

                        


#JL.model <- lm(AT ~ pop+pop:S_38836669, data = geno.and.pheno) 
JL.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno)    # linear model of trait modeled by pop, pop:SNPs specifying data in geno.and.pheno

pop.by.marker.effects <- summary(JL.model)$coefficients[-c(1:25),]      # obtain descriptive statistics for linear model and model coefficients for pop and SNP within pop, access individual columns with model coefficients
         #changed 26 back to 25 - THIS CONTRIBUTED TO THE LOSS OF POP1:snp1 EFFECT IN ALL POP.BY.MKR.EFFECT FILES
head(pop.by.marker.effects)                                              #header for this data matrix pop term, estimate, stderr, tvalue, prob, FDR adj pvalue

#Perform the FDR procedure on the P-values from pop.by.marker.effects
 res <- mt.rawp2adjp(pop.by.marker.effects[,4], "BH")
 adjp <- res$adjp[order(res$index), ]
                                       
 pop.by.marker.effects <- cbind(pop.by.marker.effects, adjp[,2])
 colnames(pop.by.marker.effects)[5] <- "FDR_Adjusted_P-values"

#Write out the results, put it in the same directory as the JL results from TASSEL
 pop.by.marker.effects <- cbind(rownames(pop.by.marker.effects), pop.by.marker.effects)
 colnames(pop.by.marker.effects)[1] <- "Term"
 
 setwd(paste(home.dir, "\\Data_for_1pct_Support_Intervals\\JL_output\\JL_model_modified_SI01\\Pop_by_mkr_est_SI01_verified\\", sep = ""))
 write.table(pop.by.marker.effects, paste("Pop.by.Marker.Effect.Estimates.from.R.", trait ,".SI01.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

#################################################################################################################################################################################
#For loop through all chromosomes, which will fit a separate GLM model for each chromosome and output all residuals

setwd(paste(home.dir,"\\Residuals_for_GWAS_SI01" ,sep = ""))
                                                                 
#test:---model.order <- model.order[-c(1:10), ]

#strategy: for loop through each chromosome. Use the "which" statement to get rid of SNPs on the chromsome being tested. The model.order object will be utilized for this analysis step.
 # Thn fit a model based off of the remaining SNPs in model.order
geno.and.pheno.mod <- geno.and.pheno[complete.cases(geno.and.pheno[,3]),]   #CBK ADDED to obtain geno.and.pheno set without missing trait BLUEs

#resid <- as.data.frame(cbind(seq(1:nrow(geno.and.pheno)),as.character(geno.and.pheno[,1])))         # has only 2 columns - order and geno_code [3432,2]
 
 #using only genos with trait BLUES
resid <- as.data.frame(cbind(seq(1:nrow(geno.and.pheno.mod)),as.character(geno.and.pheno.mod[,1])))

for(j in 1:10){
  model.order.temp <- model.order
  if(length(which(model.order.temp[,2] == j))>0){
   model.order.temp <- model.order.temp[-which(model.order.temp[,2] == j),]
  } #end if
  index <- model.order.temp[,5]+3
  base.model <- paste(trait, " ~ pop+",sep = "")
  for(k in index){
   base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno.mod)[k],sep = "")
  } #end for(k in 4:ncol(geno.and.pheno))

  #Fit the model specific to the chromosome being tested                                                     
  Residual.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno.mod)    #dim(geno.and.pheno... {3432,16}
  #colnames(resid) <- c("V3","V2","V1")             #added by CBK
  
  resid.chr <-as.data.frame(as.matrix(resid(Residual.model)))     #only has V1, residuals for linear model of trait modeled by marker by pop [3314,1]
  #resid.chr <-cbind(geno.and.pheno.mod, as.data.frame(as.matrix(resid(Residual.model) ) ) )

  #resid <- merge(resid, resid.chr, by.x = "V1", by.y = "row.names")
  resid <- cbind(resid, resid.chr)
 

}#end for(j in 1:10)

  resid <- resid[,-1]
  colnames(resid) <- c("Sample", "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10")

 write.table(resid, paste("Resid.", trait ,".SI.01.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

}# End for(i in trait)




