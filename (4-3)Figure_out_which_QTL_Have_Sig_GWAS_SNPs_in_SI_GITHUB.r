rm(list = ls())
##########################################################################

###Required files:
### (1) Tabular summary file with overlapping support intervals generated in script (4-2)
### (2) candidate gene list (see sample)
### (3) GWAS Results Summaries generated from script (3-1), "Complete Results..." file

extractinate.GWAS.SNPs <- function(){
  #Function 1 should begin here
  
  #Obtain chromosome, start and stop support interval positions, from the "tabular summary"
  row.of.interest <- which((tabular.summary[,1] == trait)  & (tabular.summary[,6] == QTL))
  chr <- tabular.summary[row.of.interest, 2]
  QTL.start.bp <- tabular.summary[row.of.interest, 4]
  QTL.stop.bp <- tabular.summary[row.of.interest, 5]
  
  
  #To be rigorous, add an object that has all of the names of the NAM founders
  pop.seq <- as.data.frame(as.factor(c("pop01", "pop02", "pop03", "pop04", "pop05", "pop06", "pop07", "pop08", 
                                       "pop09", "pop10", "pop11", "pop12", "pop13",
                                       "pop14", "pop15", "pop16", "pop18", "pop19", 
                                       "pop20", "pop21", "pop22", "pop23", "pop24", 
                                       "pop25", "pop26")))
  founder.names <- as.data.frame(c("B97", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", 
                                   "CML69", "HP301", "Il14H", "KI11", "KI3", "KY21", "M162W", "M37W", "MO18W", 
                                   "MS71", "NC350", "NC358", "OH43", "OH7B", "P39", "TX303", "TZI8"))
  NAM.pops <- cbind(pop.seq, founder.names)
  colnames(NAM.pops) <- c("Pop.num", "Pop.Founders")
  
  #Obtain the GRZM IDs of genes that are in the support interval
 
   cand.chr <- candidate.gene.list[which(candidate.gene.list[,2] == chr),] 
   
   cand.gene.names <- cand.chr[ which(( (cand.chr[,3] > QTL.start.bp) & (cand.chr[,4] < QTL.stop.bp) ) |  
             ( (cand.chr[,3] < QTL.start.bp) & (cand.chr[,3] > QTL.stop.bp)  ) |
             ( (cand.chr[,4] < QTL.stop.bp) & (cand.chr[,4] > QTL.stop.bp) ) ), 2]
   
   genes.identified.in.SI <- NULL
   for(k in 1:length(cand.gene.names)) genes.identified.in.SI <- paste(genes.identified.in.SI,", ", cand.gene.names[k], sep = "")
  #Final product of this phase: genes.identified
 
  
  
  #Obtain the GWAS results for this trait, and filter out all SNPs with RMIP < cutoff specified at top of script (5 is default)
  #if (trait %in% multicollinear.traits){
  #GWAS.results <- read.table(paste(home.dir,location.of.GWAS.results,trait,"_no_collin\\Complete_Results_",trait,"_no_collin.txt", sep = ""),head = TRUE)} else {
  GWAS.results <- read.table(paste(raw.GWAS.results.path,trait,"_iterations/Complete_Results_",trait,".txt", sep = ""),head = TRUE)#}
  GWAS.results <- GWAS.results[-which(GWAS.results[,ncol(GWAS.results)] < RMIP_cutoff),]
  
  
  #If there are no SNPs within the JL support interval, STOP the script. Print out a message saying that JL and GWAS results do not overlap.
  GWAS.results.within.SI <- GWAS.results[which((GWAS.results[,1] == chr) & (GWAS.results[,2] > QTL.start.bp) & 
                                         (GWAS.results[,2] < QTL.stop.bp)),]
  
  number.of.SNPs.in.SI <- nrow(GWAS.results.within.SI)
  
  return(list(genes.identified.in.SI = genes.identified.in.SI, number.of.SNPs.in.SI = number.of.SNPs.in.SI))
}#end extractinate.GWAS.SNPs 


##########################################################################
##########################################################################
##########################################################################



#Set the working directory
setwd("C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/JL/")                                #CHD added 5-27
home.dir <- getwd()
raw.GWAS.results.path = paste(home.dir, "/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/GWAS_Analysis/GWAS_25fam_HMPonly_TASSEL3_alpha01_2015_corr/",sep='')
trait.set <- "tocos" #Options are "carots" or "tocos"
RMIP_cutoff = 5
geno.path = paste(home.dir, "/validate_CBK.AEL/",sep='')
tabSummary.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Tabular_Summaries/",sep='')

#Read in the appropriate files
# summary without common support intervals
tabular.summary <- read.table(paste(tabSummary.path,"/Tab_Sum_with_Common_SI_Info_left_bound.txt", sep = ""), head = TRUE)
# summary with common support intervals
#tabular.summary <- read.table(paste(home.dir,"\\(16)Generating Robust Files for Group Review\\Generating overlapping support intervals\\Tab_Sum_T3_2015_SI01_with_Common_SI_Info.txt", sep = ""), head = TRUE)

#Read in the candidate genes
if (trait.set == "carots"){
  #setwd(paste(home.dir, "\\(16)Generating Robust Files for Group Review\\Adding family info to JL tabulated results_SI01", sep = ""))
  candidate.gene.list <- read.table(paste(geno.path,"carot_candidate_genes.txt",sep=''), head = TRUE)
}
if (trait.set == "tocos"){
  candidate.gene.list <- read.table(paste(geno.path,"Tocochromanol_Candidate_Gene_List_GRZMs_R.formatted.txt",sep=''), head = TRUE)
}

#multicollinear.traits <- c("BCAR_RUV", "ZEA_RUV", "TOTCAR_RUV")

#Run this for loop
the.GRZM.results <- NULL
the.SNP.results <- NULL
for(i in 1:nrow(tabular.summary)){
  trait <- as.character(tabular.summary[i,1])
  if(trait == "Total_Tocopherols"){trait = "totalT"}
  if(trait == "Total_Tocotrienols"){trait = "totalT3"}
  if(trait == "Total_Tocochromanols"){trait = "totalTocochrs"}
  QTL <- as.character(tabular.summary[i,6])
  #For loop through all rows of the tabluar summary
  the.GWAS.SNP.numbers.and.GRZMs.in.SI <- extractinate.GWAS.SNPs()
  the.GRZM.results <- c(the.GRZM.results, the.GWAS.SNP.numbers.and.GRZMs.in.SI$genes.identified.in.SI)
  the.SNP.results <- c(the.SNP.results, the.GWAS.SNP.numbers.and.GRZMs.in.SI$number.of.SNPs.in.SI)
}
######################################################################
#Append the.GRZM.results and the.SNP.results to the tabular.summary
tabular.summary <- cbind(tabular.summary, the.GRZM.results, the.SNP.results)

#Output the updated tabular summary
write.table(tabular.summary, paste(tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_SI_with_Overlapping_GWAS_SNPs_2015.txt",sep=''), sep = "\t", row.names = FALSE, quote = FALSE)
