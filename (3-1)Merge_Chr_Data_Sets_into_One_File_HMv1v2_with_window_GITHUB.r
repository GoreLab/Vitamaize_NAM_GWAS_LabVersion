#Remove all objects (i.e., data sets) from R
rm(list=ls())

###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) GWAS model files for each chromosome and bootstrap sample



#For creating windows based on genetics map: Read in the SNP Name (containing physical cooordinates in its name), Genetic postion (cM), and Chromosome
genotypes <- read.table("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(9)JL Analysis\\Permutations\\GBS_SNPs\\imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt", head = TRUE)
#Extract only the first five columns of genotypes for creating windows
window.data <- genotypes[,1:5]


#Set our working directory


setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(10)GWAS Analysis\\RUV GWAS 25fam_HMPonly_TASSEL3_alpha01_2015 corr\\")

working.directory <- "C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(10)GWAS Analysis\\RUV GWAS 25fam_HMPonly_TASSEL3_alpha01_2015 corr\\"


#Input a list of paths, and use the "List" function to combine them into one object

#trait<- c("LUT_RUV","ZEA_RUV","ZEI_RUV","BCRY_RUV",
#"ACAR_RUV","BCAR_RUV","PHYF_RUV","THLYC_RUV",
#"LUT_MSS","ZEA_MSS","ZEI_MSS","BCRY_MSS",
#"ACAR_MSS","BCAR_MSS","PHYF_MSS","THLYC_MSS",
#"TOTCAR_RUV","TOTCAR_MSS") 
          window.size.genetic <- 1 #Give the Window size in cM
num.chr <- 10
trait <- c("ACAR_RUV", "BCRY_RUV", "PHYF_RUV", "THLYC_RUV", "ZEI_RUV")
#i <- trait

#For each path (i.e., for each trait)

for(i in trait){
condensed.data <-NULL

  path <- working.directory
  
  #setwd(paste(working.directory, subdirectory.vector[[i]], sep = ""))

#Read in the ten files, and summarize their counts

     for(j in 1:10){
   #for(j in 1:num.chr){
      #Read in the data
      iter.seq <- c(0,10,20,30,40,50,60,70,80,90)
      data.chr <- NULL
      for(k in iter.seq){
        #data.path <- paste(path,i,"\\",i,"_model_chr",j,"_iter",k,".txt", sep = "")
        #data.temp <- try(read.table(data.path, head=FALSE, sep = ""))
        data.temp <- try(read.table(paste(path,i,"\\",i,"_model_chr",j,"_iter",k,".txt", sep = ""), head = FALSE))
        if(!inherits(data.temp, "try-error")) data.chr <- rbind(data.chr, data.temp)
      } #end  for(k in iter.seq)
      #Keep the chromosome and bp position
      data.chr <- data.chr[,c(1:2,4)]
            
      #"D:\\NAM_Carot_GWAS_v1v2\\TRAIT\\TRAIT_model_chr1_iter"
    
      #Sort the data by Chromosome position
              
      if(length(data.chr)> 0) data.chr <- data.chr[order(data.chr[,2]),]
     
     
      # To ensure that a locus that made it into the final model when coded as 1 = missing, 0 = otherwise is not being counted also when it is in the model as coded as a 1= minor, 0 = missing,
      # the following loop splits the counting by "allele" types
 
      for(k in 1:length(unique(data.chr[,2]))){
        data.chr.tmp <- data.chr[which(data.chr[,2]==unique(data.chr[,2])[k]  ),]
        #Calculate the RMIPs
        # The [,-4]  here is to get rid of the fourth column,, which is a duplicate of the second column (bp position). which gets automatically added from the table() function
        RMIP.times.100 <-  cbind(rep(j, length(unique(data.chr.tmp[,2]))), unique(data.chr.tmp[,c(2,3)]), table(data.chr.tmp[,2]))[,-4] 
        #Append the results to a data set for all chromosomes
        condensed.data <- rbind(condensed.data, RMIP.times.100)
        rm(RMIP.times.100)
        rm(data.chr.tmp)
        
      } #End for(k in 1:length(unique(data.chr[,3])))
   
   
      #Obtain Frequencies
   
       
    
    rm(data.chr)
   
   } #end for(j in 1:10)
   
   
    colnames(condensed.data) <- c("Chr", "bp", "allele", "RMIPx100" )


#Sort the data from largest to smallest RMIP
    condensed.data <- condensed.data[order(condensed.data[,4], decreasing = TRUE),]

#Output them into a single file                                     

   write.table(condensed.data, paste(path,i, "\\Complete_Results_", i,".txt", sep = "" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



######################################################
#Now, sum up all RMIP scores in evenly distributed windows across each chromosome
results.within.windows <- NULL
for(j in 1:num.chr){
  condensed.within.chr <- condensed.data[which(condensed.data[,1] == j),]
  # Get all of the window data of the ith chromosome
  window.data.temp <- window.data[which(window.data[,3] == j),] 

  #Obtain the markers that are only at the prespecified "window.size.genetic" cM intervals
  window.intervals <-round(c(seq( round(min(window.data.temp[,5]),0), round(max(window.data.temp[,5]),0), by = window.size.genetic)),1)
  #Make sure that the extreme ends of the chromosomes are note excluded.
  if(min(window.data.temp[,5]) < round(min(window.data.temp[,5]),0)){
     window.intervals <- c(min(window.data.temp[,5]), window.intervals)} else{
     window.intervals[1] <- min(window.data.temp[,5])
  }   
  if(max(window.data.temp[,5]) > round(max(window.data.temp[,5]),0)){
     window.intervals <- c(window.intervals, max(window.data.temp[,5]))} else{
     window.intervals[length(window.intervals)] <- max(window.data.temp[,5])
  }
 
 
  #Obtain the physical coordinates of these windows
  window.intervals.physical <- window.data.temp[which(round(window.data.temp[,5],1)%in%window.intervals), 4]
 
  if(which(round(window.data.temp[,5],1)%in%window.intervals)[1] != 1) window.intervals.physical <- c(window.data.temp[1, 4], window.intervals.physical)
 
  if(which(round(window.data.temp[,5],1)%in%window.intervals)[length(which(round(window.data.temp[,5],1)%in%window.intervals))] != nrow(window.data.temp)) window.intervals.physical <- c(window.intervals.physical, window.data.temp[nrow(window.data.temp), 5])
 
  if(length(window.intervals.physical) != length(window.intervals)){
    print(paste("There was an error with obtaining intervals of length ", window.size.genetic, " cM on Chromosome ", 2, sep = ""))
    break
  }
  
  results.within.windows.chr <- NULL
  for(k in 1:(length(window.intervals.physical)-1)){
    data.int <- condensed.within.chr[which( (condensed.within.chr[,2] > window.intervals.physical[k] )  & (condensed.within.chr[,2] < window.intervals.physical[k+1] )  ),]
    if(nrow(data.int) !=0 ) results.within.windows.chr <- rbind(results.within.windows.chr, c(j, ((window.intervals.physical[k+1] + window.intervals.physical[k])/2),  sum(data.int[,4]) )  )  
  }

  results.within.windows <- rbind(results.within.windows, results.within.windows.chr)
 

} # End    for(j in 1:num.chr)


# Go to your permutation test results to determine the denomiator degrees of freedom

    colnames(results.within.windows) <- c("Chr", paste("Midpoint_in_",window.size.genetic ,"_cM_window", sep = ""), "Sum" )


#Sort the data from largest to smallest RMIP
    results.within.windows <- results.within.windows[order(results.within.windows[,3], decreasing = TRUE),]

#Output them into a single file                                     

   write.table(results.within.windows, paste(path,i, "\\Complete_Results_within_ME_Windows_", i,".txt", sep = "" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

 rm(results.within.windows)
 rm(condensed.data)
 
}#end for(i in trait)                                                         


#############################################
# FYI:  This error message is fine... 
#simply indicates that there are sampling 
#iterations of particular chromosomes 
#which do not contain significant associations
#Error in read.table(paste(path, i, "\\", i, "_model_chr", j, "_iter",  : 
  #no lines available in input
  #
  #Script has been adjusted by AEL on Nov 6 2013 to cycle through the loop even if there are files with zero data in them
  ##########################################




#  Obtain the P-value from an F distribution
pf(Test.Statistic, df1 = 26, df2 = 4552 )
