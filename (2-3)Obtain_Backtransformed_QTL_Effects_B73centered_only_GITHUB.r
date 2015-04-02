 ###Required files:
### (1) lambda values from Box-Cox procedure for all traits
### (2) broad-sense heritabilities, line mean basis
### (3) Marker effect estimates from script (2-2)
### (4) untransformed studentized BLUEs


back.transform.est <- function(data, data.lambda, trait, est.column = 2){

 Back.Trans.Est <- NULL
 for(i in 1:nrow(data)){
  constant <- 1
  #Obtain the alleleic effect estimates
   est.trans <- data[i,est.column]
 
  #Obtain the the appropriate lambda
  # original lambda <- data.lambda[which(as.character(data.lambda[,1])==trait),2] #Begin Here
   lambda <- data.lambda[which(as.character(data.lambda[,2])==trait),3] #Begin Here
 
  #Obtain 1 Plus effect Estimate
  trans.plus.1 <- constant+est.trans
  #If trans.plus.1 is < 0 and lambda is < 0, then this function will return an NA on line 26. Thus, we have the while loop before
  while(trans.plus.1 <= 0){
    constant = constant+1
    trans.plus.1 <- constant+est.trans
  }#end while
  
  
  #Backtransform. Have logical statements to separate when lambda is equal to zero
  if(lambda!=0){
   back.trans.1 <- constant^(1/lambda)
   est.backtrans <- trans.plus.1^(1/lambda)
  }else{
   back.trans.1 <- exp(constant)
   est.backtrans <- exp(trans.plus.1)  
  }
 

 #Obtain the backtransformed effect estimate, and append it to the series of results
 effect.backtrans <- est.backtrans - back.trans.1
 
 #Append the results
 Back.Trans.Est <- c(Back.Trans.Est, effect.backtrans)
 } # end for(i in 1:nrow(data))
 return(list(lambda = lambda, Back.Trans.Est=Back.Trans.Est))

}#end back.transform.est


####################################################################################


#Set the working directory
setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
home.dir <- getwd()
BLUE.or.BLUP <- "BLUE"


#Set the output directory
output.directory <- "Backtransformed"


#list all of the traits
traits <- c("LUT_RUV", "ZEI_RUV", "BCRY_RUV", "ACAR_RUV",  "PHYF_RUV", "THLYC_RUV") 
#trait.multicollinearity <- c("BCAR_RUV","TOTCAR_RUV")
#traits <- trait.multicollinearity

  #note that for multicollinear traits, the pop by marker estimate files will contain .no.multicollinearity.test.XXX preceding trait name

#Read in the lambda values
data.lambda <- read.table(paste(home.dir, "\\(8b)BLUE Outlier Rem and Trans Redo\\Carot_BLUE_lambda_values_all.txt", sep = ""), head = TRUE)



#For loop through the traits
for(i in traits){
  #Read in the trait data
  setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\", sep = ""))
  this.data <- read.table( paste("Pop.by.Marker.Effect.Estimates.from.R.", i,".SI01_2015.txt", sep = ""), sep = "\t", head = TRUE)
  
  
  
  #setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_1pct_Support_Intervals\\JL_output\\JL_model_modified_SI01\\Pop_by_mkr_est_SI01_verified\\", sep = ""))
  #this.data <- read.table( paste("Pop.by.Marker.Effect.Estimates.from.R.", i,".SI01.txt", sep = ""), sep = "\t", head = TRUE)
  #this.data <- read.table( paste("Pop.by.Marker.Effect.Estimates.from.R.no.multicollinearity.test.", i,".txt", sep = ""), sep = "\t", head = TRUE)
  
  
  #Back transform the estimate
  this.backtrans <- back.transform.est(data = this.data, data.lambda, i, est.column = 2) 
  
  #Append it to the trait data
  this.data <- cbind(this.data, this.backtrans$Back.Trans.Est)
  
  #Change the column name of the backtransformed column to make it more understandable
  colnames(this.data)[ncol(this.data)] <- "Back_Transformed_Alleleic_Effect_Estimate"
 
  #Export the results 
  setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\", sep = ""))
   write.table(this.data, paste("NonStd.Pop.by.Marker.Effect.Est.from.R.B73centered.SI01.2015.",i,".txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    #write.table(this.data, paste("NonStd.Pop.by.Marker.Effect.Est.from.R.B73centered.",i,".no.col.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

  
  
}#for(i in traits)




























#Old code that can be discarded

####################################################################################
#rm(list = ls())

#save.image("C:\\Users\\alex\\Desktop\\JL_Workspace.R") #Alex, save the image on your laptop at home. You know those slow VPN connections
#Use the multtest library to obtain FDR-adjusted P-values
library(multtest)

setwd("G:\\Lipka_Hal\\Tocos_NAM_2009_2010\\Joint_Linkage_Analysis")
home.dir <- getwd()
BLUE.or.BLUP <- "BLUE"

#Read in the broad-sense heritabilities on a line-mean basis
heritabilities <- read.table("G:\\Lipka_Hal\\Tocos_NAM_2009_2010\\BLUPs\\Heritabilities_Results\\Heritabilities_and_Standard_Errors_NAM_Tocos_2009_2010.txt", sep = "\t", head = TRUE)


#Read in the lambda values


#Eventually, loop around the traits. This loop will begin here

trait <- c("DT3", "GT3", "AT3", "DT", "GT", "AT") 

#trait <- "GT"

for(i in trait){


  ##############################################################################################################################################
  #Read in the QTL effects
  setwd(paste(home.dir, "\\JL_Results_", BLUE.or.BLUP,"\\Pop_by_Marker_Terms_from_R\\",  sep = ""))
  data <- read.table( paste("Pop.by.Marker.Effect.Estimates.from.R.", i,".txt", sep = ""), sep = "\t", head = TRUE)


  #Parse out the columns so that you have one column for the pops and one column for the markers
  
  popID <- substr(data[,1], 1,5)
  markerID <-  substr(data[,1], 7,10000)
  
  #Initialize a matrix of output results
  result.vector <- NULL


  #Read in the BLUEs of the 282
   BLUEs <- read.table(paste("G:\\Lipka_Hal\\Tocos_NAM_2009_2010\\BLUPs\\Formatted_BLUEs_from_ASREML\\BLUPs_",i,"_BLUE.txt", sep = ""), sep = "\t", head = TRUE)
   BLUEs.282 <- BLUEs[which(BLUEs[,1] == 27),]
   SD.BLUEs.282 <- sd(BLUEs.282[,3])
   
  #Read in the broad-sense heritability
   herit.trait <- heritabilities[which(heritabilities[,1] == i),4] 
   
   scaling.constant <- herit.trait/SD.BLUEs.282 
    
  #for each marker I am noly looking at the first 13 QTL becuase this is the minimum number of QTL calculated among the first six tocochromanol compounds
  for(j in unique(markerID)[1:13]){

  
    data.reduced <- data[which(markerID == j),]  
    popID.reduced <- popID[which(markerID == j)]
  
    #Calculate the back-transformed alleleic effect estimates and append it to the data.reduced frame
    
    
   back.trans.results <-  back.transform.est(data = data.reduced, data.lambda = data.lambda, trait = i, est.column = 2)
    
   data.reduced <- cbind(data.reduced, back.trans.results$Back.Trans.Est)
    
    #for (family k in 1:25)
    for(k in 1:25){
      #for (family l in (k+1):2 6)
      for(l in (k+1):26){
        if((k != 17) & (l != 17)){
          if(k < 10){
            a <- scaling.constant*data.reduced[which(popID.reduced == paste("pop0", k, sep = "")),ncol(data.reduced)]
            #se.a <- scaling.constant*data.reduced[which(popID.reduced == paste("pop0", k, sep = "")),3]
          }else{
            a <- scaling.constant*data.reduced[which(popID.reduced == paste("pop", k, sep = "")),ncol(data.reduced)]
            #se.a <- scaling.constant*data.reduced[which(popID.reduced == paste("pop", k, sep = "")),3]
          }
          if(l < 10){
            b <- scaling.constant*data.reduced[which(popID.reduced == paste("pop0", l, sep = "")),ncol(data.reduced)]
            #se.b <- scaling.constant*data.reduced[which(popID.reduced == paste("pop0", l, sep = "")),3]
          }else{
            b <- scaling.constant*data.reduced[which(popID.reduced == paste("pop", l, sep = "")),ncol(data.reduced)]
            #se.b <- scaling.constant*data.reduced[which(popID.reduced == paste("pop", l, sep = "")),3]
          }
          #Subtract the effect sizes
          effect.size <- a - b
          
          #Calculate the standard error
          #se.effect.size <- sqrt((se.a^2) + (se.b^2))
          
          #Append results (marker, pop_k, pop_l, effect size, standard error assuming independence)
          if((k==1)&(l == (k+1))) result.vector <- rbind(result.vector, c(j, "B73", k, round(a,4)))
          if( l == (k+1) ) result.vector <- rbind(result.vector, c(j, "B73", l, round(b,4)))
          result.vector <- rbind(result.vector, c(j, k, l, round(effect.size,4)))
        } #end  if((k != 17) | (l != 17))
        if((k==17)&(l==18)){
          b <- scaling.constant*data.reduced[which(popID.reduced == paste("pop", l, sep = "")),ncol(data.reduced)]
          #se.b <- scaling.constant*data.reduced[which(popID.reduced == paste("pop", l, sep = "")),3]
          result.vector <- rbind(result.vector, c(j, "B73", l, round(b,4)))
          print("The loop went into if((k==17)&(l==18))")
        }#end  if((k==17)&(l==18))
      }#end for(l in (k+1):26)
    }#end for(k in 1:25)
    
  #Put headers on the resutls
 
  }  #end  for(j in unique(markerID))
  
  colnames(result.vector) <- c("Marker", "Parent_1", "Parent_2", "Effect_Size")  
  
  #Export the results
  setwd(paste(home.dir, "\\Standardize_QTL_Effects\\Scaled_Results_Backtransformed_test",  sep = ""))
  write.table(result.vector, paste("Standaridized_QTL_Effects_for.", i ,".txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)



} #End  for(i in trait)









































#################################################################################################################################################################################
#The code below is old and can be deleted once this sciprt is written

#################################################################################################################################################################################
#For loop through all chromosomes, which will fit a separate GLM model for each chromosome and output all residuals

setwd(paste(home.dir,"\\Residuals_for_GWAS" ,sep = ""))

#test:---model.order <- model.order[-c(1:10), ]

#strategy: for loop through each chromosome. Use the "which" statement to get rid of SNPs on the chromsome being tested. The model.order object will be utilized for this analysis step.
 # Thn fit a model based off of the remaining SNPs in model.order
resid <- as.data.frame(cbind(seq(1:nrow(geno.and.pheno)),as.character(geno.and.pheno[,1])))

for(j in 1:10){
  model.order.temp <- model.order
  if(length(which(model.order.temp[,2] == j))>0){
   model.order.temp <- model.order.temp[-which(model.order.temp[,2] == j),]
  } #end if

  if(nrow(model.order.temp) != 0){
    index <- model.order.temp[,5]+3
    base.model <- paste(i, " ~ pop+",sep = "")
    for(k in index){
      base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")
    } #end for(k in 4:ncol(geno.and.pheno))
  } else{
    base.model <- paste(i, " ~ pop",sep = "")
 
  } 
 
  #Fit the model specific to the chromosome being tested
  Residual.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno) 

  
  resid.chr <-as.data.frame(as.matrix(resid(Residual.model)))

  resid <- merge(resid, resid.chr, by.x = "V1", by.y = "row.names")

}#end for(j in 1:10)

  resid <- resid[,-1]
  colnames(resid) <- c("Sample", "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10")

 write.table(resid, paste("Resid.", i ,".txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)



}# End for(i in trait)




