rm(list = ls())

#Use the multtest library to obtain FDR-adjusted P-values

###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) candidate gene list (see sample)
### (3) transformed BLUEs
### (4) JL model output, with NA in pop term row (modified version)
### (5) Marker effect estimates from script (2-2)


setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
home.dir <- getwd()

#Read in the gneotypic data
setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\GBS_SNPs\\", sep = ""))
genotypes <- read.table("imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt", head = TRUE)

#Read in the candidate genes
setwd(paste(home.dir, "\\(16)Generating Robust Files for Group Review\\Adding family info to JL tabulated results_SI01", sep = ""))
candidate.gene.list <- read.table("carot_candidate_genes.txt", head = TRUE)

# NOTE that trait loop was not run... rather each trait was run individually

#Loop around the traits.
trait <- c("LUT_RUV", "ZEA_RUV","ZEI_RUV","BCRY_RUV","ACAR_RUV",
  "PHYF_RUV","THLYC_RUV", "BCAR_RUV", "TOTCAR_RUV") 
#trait.multicollinearity <- c("BCAR_RUV", "TOTCAR_RUV",  "ZEA_RUV") 
BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"
tocos.or.carots <- "Carots" #Options are "Tocos" or "Carots"
test <- FALSE #Options are TRUE and FALSE
                                
pop.seq <- as.data.frame(as.factor(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13",
                                    "14", "15", "16", "18", "19", "20", "21", "22", "23", "24", "25", "26")))
founder.names <- as.data.frame(c("B97", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", 
                                 "CML69", "Hp301", "Il14H", "Ki11", "Ki3", "Ky21", "M162W", "M37W", "Mo18W", 
                                 "MS71", "NC350", "NC358", "Oh43", "Oh7B", "P39", "Tx303", "Tzi8"))
NAM.pops <- cbind(pop.seq, founder.names)
colnames(NAM.pops) <- c("Pop.num", "Pop.Founders")
#trait <- trait.multicollinearity

Results <- NULL
for(i in trait){

  #Read in the results from JL analysis
  
  #if(i %in% trait.multicollinearity){
    #setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\JL_output\\JL_model_modified_SI01\\",  sep = ""))
    #TASSEL.model.results <- read.table(paste(i, "_JL_model_modified_SI01_2015.txt", sep = "") , head = TRUE)  
  #}else{
    setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\JL_output\\JL_model_modified_SI01\\",  sep = ""))
    TASSEL.model.results <- read.table(paste(i, "_JL_model_modified_SI01_2015.txt"   , sep = "") , head = TRUE)
  #}
  
  if(is.numeric(TASSEL.model.results[,2]) == FALSE){
    print(paste("The chromosome column in ", i," model results needs to be numeric. Please change this before proceeding.", sep = ""))
    break;
  }


  #Obtain corresponding physical bp positions of the support intervals
  #I am going to use a for loop so that all SNPs on the same chromosome are isolated. This will prevent selecting a SNP on a different chromosome with the same genetic position.
  support.int.physical <- NULL
  for(j in 2:nrow(TASSEL.model.results)){
   geno.chr <- genotypes[which(genotypes[,3] %in% TASSEL.model.results[j,2]),1:5] 
   support.int.left.pos <- geno.chr[which(round(geno.chr[,5],1) == round(TASSEL.model.results[j,11],1)), 4]
   support.int.right.pos <- geno.chr[which(round(geno.chr[,5],1) == round(TASSEL.model.results[j,12],1)), 4]
   support.int.physical <- rbind(support.int.physical, c(support.int.left.pos, support.int.right.pos))
  }
    
  #Determine if any of the 60 candidate genes are in the support intervals
  cand.gene <- rep(NA, nrow(support.int.physical))

  for(j in 1:nrow(support.int.physical)){
   cand.chr <- candidate.gene.list[which(candidate.gene.list[,2] == TASSEL.model.results[(j+1),2]),] 
   
   cand.gene.names <- cand.chr[ which(( (cand.chr[,3] > (support.int.physical[j,1])) & (cand.chr[,4] < (support.int.physical[j,2])) ) |  
             ( (cand.chr[,3] < (support.int.physical[j,1])) & (cand.chr[,3] > (support.int.physical[j,2]))  ) |
             ( (cand.chr[,4] < (support.int.physical[j,2])) & (cand.chr[,4] > (support.int.physical[j,2])) ) ), 1]
   
   genes.identified <- NULL
   for(k in 1:length(cand.gene.names)) genes.identified <- paste(genes.identified,", ", cand.gene.names[k], sep = "")
   
   cand.gene[j] <- genes.identified 
  }
  
  
    
  ########################## Calculate PVE
  #### Read in the phenotypic data, and merge it to the genotypic data
  setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Phenotypes\\Phen_files_no_trait_prefix\\",  sep = ""))
  pheno.data <- read.table(paste(i,"_BLUE_outtrans.txt", sep = ""), head = TRUE)
  
  geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)]
  geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
  colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))



  #pheno.data will always have more data becuase IBM is included in the phenotypic data.

  geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

  #Add a population column
  geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
  colnames(geno.and.pheno)[2] <- "pop"

  
  ####Extract the trait value, RIL ID, and family
  pheno.data.temp <- geno.and.pheno[,1:3]
  
  
  #Remove all missing phenotypic data
  
  if(length(which(is.na(pheno.data.temp[,3])))>0) pheno.data.temp <- pheno.data.temp[-which(is.na(pheno.data.temp[,3])),]
  
  #For validation. 
  if(test == TRUE){
    setwd(paste(home.dir, "\\PVE_Documentation_and_Validation\\",  sep = ""))
    write.table(pheno.data.temp, paste("Pheno.data.temp.", i,".txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  ####Calculate the following summary statistics: 
  #Variance within each family

  Vf.vector <- as.vector(tapply(pheno.data.temp[,3], pheno.data.temp[,2], var))

  #Sample size within each family
  nf.vector <- as.vector(tapply(pheno.data.temp[,3], pheno.data.temp[,2], length))
 
  #Overall sample size
  N <- sum(nf.vector)
  nf.vector.by.N <- nf.vector/N
 
  den <- (1/N)*crossprod(nf.vector, Vf.vector)

    
  #Read in the marker data
  #if(i %in% trait.multicollinearity){
    #setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_1pct_Support_Intervals\\JL_output\\JL_model_modified_SI01\\Pop_by_mkr_est\\",  sep = ""))
    #pop.by.marker.effects <-  read.table(paste("Pop.by.Marker.Effect.Estimates.from.R.no.multicollinearity.test.", i,".txt", sep = ""), sep = "\t", head = TRUE)
  #}else{
    #setwd(paste(home.dir, "\\((9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\",  sep = ""))
    setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\")
    
    pop.by.marker.effects <-  read.table(paste("Pop.by.Marker.Effect.Estimates.from.R.", i,".SI01_2015.txt", sep = ""), sep = "\t", head = TRUE)
   #pop.by.marker.effects <-  read.table(paste("Pop.by.Marker.Effect.Estimates.from.R.",i,".SI01_2015.txt", sep = ""), sep = "\t", head = TRUE)
  #}

  
  #Separate out the pop and maker terms and append them to "pop.by.marker.effects"
  pop.term <- as.factor(substr(as.character(pop.by.marker.effects[,1]), start = 4, stop = 5))
  marker.term <-  as.factor(substr(as.character(pop.by.marker.effects[,1]), start = 7, stop = nchar(as.character(pop.by.marker.effects[,1]))))
  
  pop.by.marker.effects <- cbind(pop.term, marker.term, pop.by.marker.effects)
  
  #This following sequence is used to figure out which family has no pop*marker effect estimate (this arises when the heritabilities are zero)
  pop.seq <- as.data.frame(as.factor(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13",
                         "14", "15", "16", "18", "19", "20", "21", "22", "23", "24", "25", "26")))
  
  colnames(pop.seq)[1] <- "Pop"
  
  #Set a PVE results vector to NULL
  PVE.results <- NULL
  #For loop using j as an index
  for(j in unique(pop.by.marker.effects[,2])){
    #Extract the corresponding additive effect estimates 
    additive.effects <- pop.by.marker.effects[which(pop.by.marker.effects[,2] == j), c(1,2,4)]
    
    #For validation
    if(test == TRUE){
      setwd(paste(home.dir, "\\PVE_Documentation_and_Validation\\",  sep = ""))
      write.table(additive.effects, paste("Additive.effects.temp.", i,".txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  
   
    if(nrow(additive.effects) != 25){
     #Objective: Put nf.vector, Vf.vector into a matrix, and then merge it with the additve effect estimates
     vectors <- as.data.frame(cbind(pop.seq, nf.vector, Vf.vector))    
     add.and.vectors <- merge(additive.effects, vectors, by.x = "pop.term", by.y = "Pop")
      
     #Now the nf.vectora nd Vf.vector will not have the populations with no additive effect estimates
     nf.vector.temp <- add.and.vectors[ ,4]
     Vf.vector.temp <- add.and.vectors[ ,5]    
     
     a.sq <- additive.effects[,3]^2
     num.left.term  <- ((2/N)*(crossprod(nf.vector.temp,a.sq)))
     
     num.right.term <- ((1/N)*(crossprod(nf.vector.temp,additive.effects[,3])))^2
     
     
     num <- num.left.term - num.right.term

     N.temp <- sum(nf.vector.temp)

 
     den <- (1/N.temp)*crossprod(nf.vector.temp, Vf.vector.temp)

    
     PVE.marker <- c(j, (num/den))

     
    }else{
        
    #Calculate the PVE
     a.sq <- additive.effects[,3]^2
     num.left.term  <- ((2/N)*(crossprod(nf.vector,a.sq)))
     
     num.right.term <- ((1/N)*(crossprod(nf.vector,additive.effects[,3])))^2
     
     
     num <- num.left.term - num.right.term

    
     PVE.marker <- c(j, (num/den))
    }
    #Append each result. Note that the name of the marker MUST be included as well. 
    PVE.results <- rbind(PVE.results,  PVE.marker)
       
  } #End for(j in unique(pop.by.marker.effects[,2]))
  
  
  
  ########################### End Calculate PVE 
  
  #################################Determine the number of QTL in each family, and which families have the QTL
  number.of.familes.with.QTL <- NULL
  specific.families.with.QTL <- NULL
  for(j in unique(marker.term)){
    #Obtain a subset of pop.by.marker.effects for the jth QTL only
    pop.by.marker.effects.subset1 <- pop.by.marker.effects[which(marker.term == j),]
    pop.term.subset1 <- pop.term[which(marker.term == j)]
    
    #Obtain a second subset of the pop.by.marker that are statistically significant at 5% FDR. These will be
    pop.by.marker.effects.subset2 <- pop.by.marker.effects.subset1[which(pop.by.marker.effects.subset1[,ncol(pop.by.marker.effects.subset1)] <= 0.05),]
    pop.term.subset2 <- pop.term.subset1[which(pop.by.marker.effects.subset1[,ncol(pop.by.marker.effects.subset1)] <= 0.05)]
    
    #Count the number of families that have the QTL
    number.of.familes.with.QTL <- c(number.of.familes.with.QTL, length(pop.term.subset2))
    
    #Identify the families that have the QTL
    as.data.frame(NAM.pops)
    as.data.frame(pop.term.subset2)
    families.with.QTL <- merge(as.data.frame(pop.term.subset2), as.data.frame(NAM.pops), by.x = "pop.term.subset2",
                               by.y = "Pop.num")
    families.with.QTL.as.string <- NULL
    for(k in 1:nrow(families.with.QTL)) families.with.QTL.as.string <- paste(families.with.QTL.as.string,",",families.with.QTL[k,2])
    specific.families.with.QTL <- c(specific.families.with.QTL, families.with.QTL.as.string)
    
  }#end for(i in unique(marker.term))
  
  
  
  
  
  #########################End of "Determine the number of QTL in Each family"
  
  #Make sure you have the following information: Name of Marker, peak physical position, support interval, P-value, PVE, Candidate genes in suppoprt interval

  results.trait <- cbind(TASSEL.model.results[-1,c(1,2)],  substr(as.character(TASSEL.model.results[-1,4]), start = 3, stop = nchar(as.character(TASSEL.model.results[-1,4]))), support.int.physical,
                         TASSEL.model.results[-1,c(4,9)], PVE.results[,2],cand.gene, number.of.familes.with.QTL, specific.families.with.QTL) 
      

  #Append these results to the results file
  Results <- rbind(Results, results.trait)

}# End for(i in trait)



setwd(home.dir)                                                                                                                                                                              

#Add headers to the Results
colnames(Results) <- c("Trait", "Chr", "Peak_Position", "Supp_Int_Left", "Supp_Int_Right", "Peak_Marker_Name", "P-value", "PVE", "Candidate_Gene_in_Support_Interval", "Number_of_Families_with_QTL", "Families_with_QTL")


#Write the results into a table
setwd(paste(home.dir,"\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\", sep = ""))    
#write.table(Results, paste("Tabular_Summary_of_JL_",tocos.or.carots,"_Results_for_all traits_SI01.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Results, paste("Tabular_Summary_of_JL_",tocos.or.carots,"_Results_for_all_traits_SI01_2015_T3.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)


