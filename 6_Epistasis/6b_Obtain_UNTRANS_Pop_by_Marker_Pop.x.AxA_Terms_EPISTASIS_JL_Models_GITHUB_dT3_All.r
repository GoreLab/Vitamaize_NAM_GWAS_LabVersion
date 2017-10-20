rm(list = ls())

###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version) and residual row removed from file

#Use the multtest library to obtain FDR-adjusted P-values
setwd("/Users/anybody/Documents/Toco_NAM/")
home.dir <- getwd()
geno.path = paste(home.dir,"/Geno.Pheno_Inputs/",sep='')
pheno.path = geno.path
UNTRANS.pheno.path = "/Users/anybody/Documents/Toco_NAM/Geno.Pheno_Inputs/BLUES_Untransf/"
dT3.new.path = paste(home.dir,"/Toco_NAM/Methods/dT3_removeExtremeVal_test/",sep='')
EPISTASIS.JL.path = paste(home.dir,"/Epistasis/Final_Epist_Models_forR/",sep='')
popByMarker.path = paste(home.dir, "/Epistasis/Allelic_Effect_Estimates/", sep = "")
library(multtest)

#Read in the genotypic data
genotypes <- read.table(paste(geno.path,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Eventually, loop around the traits. This loop will begin here                               
trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs") 
trait.collinear <- c("aT3","dT3","gT3","totalT3","totalTocochrs")   

for(i in trait){
#Read in the trait under study
  if(i == "dT3"){
    pheno.data = read.table(paste(UNTRANS.pheno.path,"BLUEs_Untransformed_No_Out_201304019_for_R_dT3_ExtrValRemoved.txt",sep=''),head=TRUE,stringsAsFactors=FALSE)
    pheno.data$dT3 = as.numeric(pheno.data$dT3)
  }else{
    pheno.data <- read.table(paste(UNTRANS.pheno.path,"BLUEs_Untransformed_No_Out_201304019_for_R_",i,".txt", sep = ""), head = TRUE)
  }

head(pheno.data)

#Read in the results from JL analysis
TASSEL.model.results <- read.table(paste(EPISTASIS.JL.path,"Add.and_epi.JL.model.plus.PVE.at.alpha.0.01for.",i,".NAM_chd_PVEbyFam.txt", sep = "") , head = TRUE)       

#### #########################################################################################################################################
#Get the appropriate SNPs and merge the phenotypic and genotypic data together.
  
#Parse out all of the SNPs that were selected in the model
#geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)] #used for additive-only models; CHD commented out, now need to account for markers in interaction terms.
markers_in_model = NULL
all_model_terms = TASSEL.model.results[,1]
for(term in all_model_terms){
  term_split = strsplit(term,":")
  for(substr in term_split){
    if(substr=="pop"){next}
    else{markers_in_model = c(markers_in_model,substr)}
  }
}

markers_in_model = unique(markers_in_model)
print(markers_in_model)

geno.reduced <- genotypes[which(genotypes[,1] %in% markers_in_model),-c(2:5)] 

geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))

#pheno.data will always have more data becuase IBM is included in the phenotypic data.
colnames(pheno.data)[1] = "Geno_Code"
geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

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

base.model <- paste(i, " ~ pop",sep = "")                    #By CBK: paste trait and pop+   - note that tilde means "modeled as"; CHD removed "+" after pop 5-14 because + already appended with each marker
for(k in index){                                                 #BY CBK: cycle through rows in index, which is the sequence ordered by chr and position
  base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")   #BY CBK:  writing out model of trait modeled by pop plus pop by SNPs in model
} #end for(k in 4:ncol(geno.and.pheno))     

#generate vector of all pairwise combinations, for j additive effects there should be !(j-1) interaction terms
all.SIGNIF.interaction.terms.vector <- NULL
for(j in 1:(length(the.markers)-1)){
  for(k in (j+1):length(the.markers)){
    single.interaction.term <- paste("pop:",the.markers[j],":",the.markers[k], sep = "")
    all.SIGNIF.interaction.terms.vector <- c(all.interaction.terms.vector, single.interaction.term)
  } #end k loop
} #end j loop

#format all interaction vector to fit into model
the.interaction.term.model.component <- NULL
combine.model.components <- NULL
for (m in all.interaction.terms.vector[1:(length(all.interaction.terms.vector)-1)]){
  the.interaction.term.model.component <- paste(m,"+")
  combine.model.components <- paste(combine.model.components, the.interaction.term.model.component)
}  # end adding interaction terms with the exception of the last term in the vector

#add last term in vector
final.interaction.term.model.component <- paste(combine.model.components, all.interaction.terms.vector[length(all.interaction.terms.vector)])

#Generate complete model with all additive and possible epistatic terms
#the.tested.model <- paste(base.model, "+",final.interaction.term.model.component ,sep = "")
the.tested.model <- base.model
#} #end generate model

print(base.model)
#JL.model <- lm(AT ~ pop+pop:S_38836669, data = geno.and.pheno) 
JL.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno)    # linear model of trait modeled by pop, pop:SNPs specifying data in geno.and.pheno

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
 
 #setwd(paste(home.dir, "\\validate_CBK.AEL\CHD_Tassel3fromSF_modified0.01\Allelic_Effect_Estimates.no.MultiColl\\", sep = ""))        #commented out CHD 5-11--moved up to top
 write.table(pop.by.marker.effects, paste(popByMarker.path,"Pop.by.Marker.Effect.Estimates.from.R.",i,".SI01_dT3.Redone_UNTRANS.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
}#end for i in trait