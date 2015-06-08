rm(list = ls())


###Required files:
### (1) transformed BLUEs
### (2) JL model output, with NA in pop term row (modified version)
### (3) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
get.geno.pheno.chr.and.cM <- function(){
  pheno.data <- read.table(paste(pheno.path,"BLUEs_No_Outliers_Transformed_all_for_TASSEL_",i,".txt", sep = ""), head = TRUE)
  
  #Read in the results from JL analysis
  #TASSEL.model.results <- read.table(paste(i, "_JL_model_SI01.txt", sep = "") , head = TRUE)
  
  #Read in the results from JL analysis
  if (trait.set == "carots"){
    #setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\JL_output\\JL_model_modified_SI01\\",  sep = ""))
    TASSEL.model.results <- read.table(paste(JL.path,i, "_JL_model_modified_SI01_2015.txt"   , sep = "") , head = TRUE)
  }
  
  if (trait.set == "tocos"){
    if (i %in% trait.collinear){   
      TASSEL.model.results <- read.table(paste(JL.path,"MC_corrected_",i,"_postRescan_R.formatted.txt", sep = "") , head = TRUE)        #CHD added loop 5-14 to handle both MC and non-MC traits
    } else {
      TASSEL.model.results <- read.table(paste(JL.path,i,"_model_3fromSF_0.01_Final_R.formatted.txt", sep = "") , head = TRUE)
    }
  }
  
  if(is.numeric(TASSEL.model.results[,2]) == FALSE){
   print(paste("The chromosome column in ", i," model results needs to be numeric. Please change this before proceeding.", sep = ""))
   break;
  }
  
  
  #Get the appropriate SNPs and merge the phenotypic and genotypic data together.
    
  #Parse out all of the SNPs that were selected in the model
  geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)]
  chr.of.JL.markers <- TASSEL.model.results[,2][-1] #These are in the order that the markers appear in the final model
  cM.of.JL.markers <- TASSEL.model.results[,3][-1] #These are in the order that the markers appear in the final model
  geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
  colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))
  
  
  
  #pheno.data will always have more data becuase IBM is included in the phenotypic data.
  colnames(pheno.data) <- c("entry", "trait")
  geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "entry", by.y = "row.names")
  
  #Add a population column
  geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
  colnames(geno.and.pheno)[2] <- "pop"
  
  #Get the model order of the model in R to be the same as that in TASSEL.
  seq <- (seq(1:nrow(TASSEL.model.results))-1)[-1]
  model.order <- cbind(seq, TASSEL.model.results[-1,c(2,3,4)])
  #Sort by chromosome and bp
  model.order <- model.order[order(model.order[,3]),]
  model.order <- model.order[order(model.order[,2]),]
  model.order <- cbind(model.order, seq(1:nrow(model.order)))
  #Sort by the first column of marker order so that the markers will be put into the model in the same order as they were selected.
  model.order <- model.order[order(model.order[,1]),]
  index <- model.order[,5]+3 #3 is added so because the first SNP is on the fourth column

  
  return(list(geno.and.pheno = geno.and.pheno, chr.of.JL.markers = chr.of.JL.markers, cM.of.JL.markers = cM.of.JL.markers,
              TASSEL.model.results = TASSEL.model.results, model.order = model.order, index = index))
}# end get.geno.pheno.chr.and.cM 

################################################
obtain.reduced.model.order.object <- function(){
    vector.of.markers.to.remove <- NULL
    for(k in 1:nrow(model.order)){ #for loop through model.order
      #If the kth marker is not on the chromosome of the tested SNP, then "next"
      if(genotypes[which(as.character(genotypes[,1]) == j),3] != model.order[k,2]){
        next
      }else{#Else
        #If the absolute value of the jth marker postion - kth marker position < window size
        if(abs(genotypes[which(as.character(genotypes[,1]) == j),5] - model.order[k,3]) < window.size.cM.on.either.side){
          vector.of.markers.to.remove <- c(vector.of.markers.to.remove, k)
        }#end if(abs(genotypes[which(as.character(genotypes[,1]) == j),5] - model.order[k,3]) < window.size.cM.on.either.side)
      }#end if(genotypes[which(as.character(genotypes[,1]) == j),3] != model.order[k,2])
    }#end for(k in 1:nrow(model.order))
    
    #Return a "reduced" model.order.object, which is the model order vector minus the rows that are within the vector
    if(length(vector.of.markers.to.remove) > 0){
      model.order.reduced <- model.order[-vector.of.markers.to.remove,]
    }else{
      model.order.reduced <- model.order
    }
    
    return(model.order.reduced)
} #end obtain.reduced.model.order.object


fit.full.model.obtain.log.L <- function(){
    #If marker j is already in the final model, do not fit it in the model again
    if(j %in% as.character(info.for.full.model$TASSEL.model.results[,4])[-1]){
      #Obtain all of the appropriate explanatory variables for the model (and in the correct order)
      #base.model <- paste(i, " ~ pop",sep = "")
      base.model <- "trait ~ pop"
      for(k in index){
        base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")
      } #end for(k in 4:ncol(geno.and.pheno))
      #Fit the full model and obtain the log likelihood function
      full.model.log.L <- logLik(lm(paste(base.model, sep = ""), data = geno.and.pheno))
      
      summary(lm(paste(base.model, sep = ""), data = geno.and.pheno))
    }else{#Else
      #Append the tested marker genomtypes to the geno.and.pheno.file
      geno.of.tested.SNPs <- genotypes[which(as.character(genotypes[,1]) == j),-c(1:5)]
      geno.of.tested.SNPs.formatted <-t(geno.of.tested.SNPs)
      #colnames(geno.of.tested.SNPs.formatted) <- as.character(t(geno.of.tested.SNPs[,1]))
      colnames(geno.of.tested.SNPs.formatted) <- j
      geno.and.pheno <- merge(geno.and.pheno, geno.of.tested.SNPs.formatted, by.x = "geno.and.pheno[, 1]", by.y = "row.names")
      #Fit the "full model":- obtain the log likelihood function value
      model.order.reduced <- obtain.reduced.model.order.object()
      index.reduced <- model.order.reduced[,5]+3
      #base.model <- paste(i, " ~ pop",sep = "") ###ADDED BY AEL ON 1/22/15 
      base.model <- "trait ~ pop"                                 
      for(k in index.reduced){
        base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")
      } #end for(k in 4:ncol(geno.and.pheno))
      #Append the tested marker to the model
      base.model <- paste(base.model,"+pop:",j,sep = "")
      #Fit the "full model": the joint linkage model plus the tested marker(family) included
      full.model.log.L <- logLik(lm(paste(base.model, sep = ""), data = geno.and.pheno))
      # as an additional explanatory variable; obtian the log likelihood function value
    }#end if(j %in% as.character(info.for.full.model$TASSEL.model.results[,4])[-1])
    return(full.model.log.L)
}#end fit.full.model.obtain.log.L

fit.reduced.model.obtain.log.L <- function(){
    #Fit the "reduced model":- obtain the log likelihood function value
    model.order.reduced <- obtain.reduced.model.order.object()
    index.reduced <- model.order.reduced[,5]+3
    #base.model.reduced <- paste(i, " ~ pop",sep = "")  
    base.model.reduced <- "trait ~ pop"    
    for(k in index.reduced){
      base.model.reduced <- paste(base.model.reduced,"+pop:",colnames(geno.and.pheno)[k],sep = "")
    } #end for(k in 4:ncol(geno.and.pheno))
    reduced.model.log.L <- logLik(lm(paste(base.model.reduced, sep = ""), data = geno.and.pheno))
    
    return(reduced.model.log.L)
}#end fit.reduced.model.obtain.log.L()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################




#Set the home directory, output directory, and results where the TASSEL alpha = 0.01 JL analysis results are stored
setwd("C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/JL/")                                #CHD added 5-14
home.dir <- getwd()
trait.set <- "tocos" #Options are "carots" or "tocos"
geno.path = paste(home.dir, "/validate_CBK.AEL/",sep='')
pheno.path = paste(home.dir, "/Phenotypes/",sep ='')
JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/R.JL.no.MultiColl/Toco_FinalModels_postMCCorrPostRescan-whenApplicable/",sep='')
tabSummary.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Tabular_Summaries/",sep='')

#Read in the markers for JL analysis
genotypes <- read.table(paste(geno.path,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Set up a vector of traits
#trait <- c("LUT_RUV", "ZEI_RUV", "BCRY_RUV", "PHYF_RUV", "THLYC_RUV") 
#trait.multicollinearity <- c("BCAR_RUV", "TOTCAR_RUV", "ZEA_RUV")
  #to change if multicollinear trait is being used
#trait <- trait.multicollinearity

#Loop around the traits.                                #CHD 5-14 loop made so that carots or tocos can be run automatically from top of script.
if (trait.set == "carots"){
  trait <- c("LUT_RUV", "ZEA_RUV","ZEI_RUV","BCRY_RUV","ACAR_RUV","PHYF_RUV","THLYC_RUV")    
  trait.collinear <- c("BCAR_RUV", "TOTCAR_RUV") 
} else if (trait.set == "tocos"){
  trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs")       #CHD added 5-11
  trait.collinear <- c("aT3","dT3","gT3","totalT3","totalTocochrs")   
}else{print("Error: no trait set selected. Please specify carots or tocos.")}

BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"
window.size.cM.on.either.side <- 3



for(i in trait){#For loop through the traits
  print(paste("------------------Starting ", i,"-------------------------------------------",sep = ""))
  #Obtain the "info.for.full.model" that includes all markers from the final JL model
  info.for.full.model <- get.geno.pheno.chr.and.cM()
  geno.and.pheno <- info.for.full.model$geno.and.pheno
  model.order <- info.for.full.model$model.order
  index <- info.for.full.model$index  
  
  count <- 1
  for(j in as.character(genotypes[,1])){#For loop through each marker
    if((j >0)&(floor(count/10)==count/10)) {print(paste("Finished ", count," out of ", nrow(genotypes)," SNPs. Trait: ", i, sep=""))}
    #Temporary code
    print(j)
    #End temporary code
    
    #Obtain the log likelihood function of the full model
    full.model.log.L <- fit.full.model.obtain.log.L()
    
    #Obtain the log likelihood function of the reduced model
    reduced.model.log.L <- fit.reduced.model.obtain.log.L()
    
    #Calculate the likelihood ratio test (LRT) statistic: D = -2ln(likelihood_reduced model) + 2ln(likelihood_full model)
    LRT <- (-2*reduced.model.log.L) + (2*full.model.log.L)
    
    #Calculate the LOD score, which is the LRT/4.61
    LOD.score <- LRT/4.61
    #Store the logL(Full), logL(Reduced), LRT statistic, and the LOD score
    the.results <- c(j, genotypes[which(as.character(genotypes[,1]) == j),3], genotypes[which(as.character(genotypes[,1]) == j),5],
                     full.model.log.L, reduced.model.log.L, LRT, LOD.score)
    if(count == 1){
      results <- the.results
    }else{
      results <- rbind(results, the.results)
    }
    count <- count+1 
    
    
    #if(count == 11){ #This code is for testing purposes
    #  print("11 runs of the loop was conducted")
    #  break
    #}
  }#end for(j in 1:nrow(genotypes))
  
  
  colnames(results) <- c("SNP_ID", "Chr", "cM", "logL_Full", "logL_Reduced", "Likelihood_Ratio_Test_Statistic", "LOD_Score")
  write.table(results, paste(tabSummary.path,"LOD_Score_Scan_",i,".txt", sep = ""),sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("------------------Finished ", i,"-------------------------------------------",sep = ""))
}#End for(i in trait)
