rm(list = ls())

###Required files:
### (1) Tabular summary file with JL/ GWAS/ common SI/ recombination freq info  (file has gone through scripts 4-1, 4-2, 4-3, 4-6)
### (2) LOD score scan files generated from script (4-4)


#Set the home directory, output directory, and results where the TASSEL alpha = 0.01 JL analysis results are stored
#setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(31) Tabular Summary Info for 2015 analysis\\")
trait.set = "tocos" #Options are 'Carot' and 'tocos'
setwd("C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/JL/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Tabular_Summaries/")
home.dir <- getwd()
JL.results.dir <- "\\Recombination_Suppressed_Regions\\"    #contains JL-GWAS commonSI recombsuppregion file
LOD.results.dir <- "/LOD_Score_Scans/" #contains LOD score data files from script (4-5)
output.dir <- home.dir  #one directory up from LOD score scans

the.JL.results <- read.table(paste(home.dir, JL.results.dir, "Tab_Sum_",trait.set,"_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions.txt", sep = ""),head = TRUE)

#initialize data vector
LOD.score.vector <- NULL

#generate loop to extract LOD score from scan file [,7], calling in file with appropriate trait name and identifying associated SNP_ID [,1] (equiv to peak marker name [,6] in the.JL.results)
#and chr [,2] (equiv to chr in [,2])
for(i in 1:nrow(the.JL.results)){
  trait <- the.JL.results[i,1]
  if(trait == "Total_Tocopherols"){trait = "totalT"}
  if(trait == "Total_Tocotrienols"){trait = "totalT3"}
  if(trait == "Total_Tocochromanols"){trait = "totalTocochrs"}
  chr <- the.JL.results[i,2]
  peak.marker.name <- as.character(the.JL.results[i,6])
  LOD.scan.file <- read.table(paste(home.dir,LOD.results.dir, "/LOD_Score_Scan_", trait, ".txt", sep = ""), head = TRUE)
  
  #retrieve LOD row of interest and score
  row.of.interest <- LOD.scan.file[which(LOD.scan.file[,1] == peak.marker.name),]
  LOD.score <- as.character(row.of.interest[7])
  
  #compile vector of LOD scores according to JL file
  LOD.score.vector <- c(LOD.score.vector, LOD.score)
} #end i loop for rows in the.JL.results

#concatenate results with original JL table
the.new.JL.results <- cbind(the.JL.results, LOD.score.vector)
write.table(the.new.JL.results, paste(output.dir,"/Tab_Sum_",trait.set,"_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions_LODscores.txt", sep = ""), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
