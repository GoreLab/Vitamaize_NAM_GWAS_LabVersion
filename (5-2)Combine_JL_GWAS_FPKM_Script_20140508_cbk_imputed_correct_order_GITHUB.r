rm(list = ls())

###Required Files:
### (1) FPKM matrix with log2 transformed data, organized by founder, named FPKM.table.by.gene.annotation.complete.founder.matrix.txt
### (2) GAPIT source files
### (3) Tabular summary file
### (4) TRANSFORMED effect estimate files, from script (2-2)
### (5) GWAS Results Summary files, from script (3-1)
### (6) QTL array file called in script (5-1)
### (7) Imputed and extracted genotype matrices from script (5-1)
### (8) Mock placeholder Hapmap genotype files for each chromosome, used as default if there are no SNPs in the selected interval


##########################################################################
###################################This function was written by Cathy Kandianis, and named by Alex Lipka
get.me.my.FPKM.values <- function(absolute.final.data.set.FPKM = NA, gene.ID = NA, gene.name = NA,
                                  print.out.results = FALSE, output.dir = NA, home.dir = NA){
  row.of.interest <- absolute.final.data.set.FPKM[(which(substr(absolute.final.data.set.FPKM[,1],1,30) == gene.ID)),]
  
  #length for each entry = 221, data starts at [,6]
  FPKM.vector <- row.of.interest[,6:ncol(row.of.interest)]
  
  founder.names <- as.matrix(c("B73", "B97", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", 
                                   "CML69", "HP301", "IL14H", "KI11", "KI3", "KY21", "M162W", "M37W", "MO17", "MO18W", 
                                   "MS71", "NC350", "NC358", "OH43", "OH7B", "P39", "TX303", "TZI8")) 
  
  as.matrix(FPKM.vector) -> FPKM.vector
  FPKM.matrix <- NULL
  test.seq <- seq(1:27)
  #for i in seq(1:length(founder.names))                                 
  for (i in test.seq){
    #after every 8 columns of data, move to a new row
    #can do this in a loop, where i is the number of iterations or founders, and multiplicative factor is 8  
    temp.line.added <- NULL
    start.of.line <- NULL
    end.of.line <- NULL
    start.of.line <- (i*8)-7
    end.of.line <- i*8
    temp.line.added <- FPKM.vector[,start.of.line:end.of.line]
    #colnames(temp.line.added) <- NULL
    FPKM.matrix <- rbind(FPKM.matrix, temp.line.added)
  }
  
  colnames(FPKM.matrix) <- c("12_DAP", "16_DAP", "20_DAP", "24_DAP", "30_DAP", "36_DAP", "root", "shoot")   
  rownames(FPKM.matrix) <- founder.names      
  
  if(print.out.results){
    setwd(paste(home.dir, output.dir, sep = ""))
    write.table(FPKM.matrix, paste("FPKM.matrix.for.",gene.name, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
  }
  return(FPKM.matrix)
}# end get.me.my.FPKM.values



##########################################################################
#Set the working directory
setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
home.dir <- getwd()
location.of.modified.GAPIT.files <- "\\(21) GAPIT source files\\"
location.of.effect.estimates <- "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\"
location.of.GWAS.results <- "\\(10)GWAS Analysis\\RUV GWAS 25fam_alldata_alpha01_2015_FINAL\\"
output.dir <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_R\\"
get.me.my.SNPs.files.dir <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_GMMS4_imputed_matrix_T3_2015\\"           # this was changed
FPKM.file.dir <- "\\(15)Correlated Expression\\"


#Source in the modified GAPIT files
setwd(paste(home.dir,location.of.modified.GAPIT.files,sep = ""))
source("GAPIT.Fragment.Modified.R")
source("GAPIT.HapMap.Modified.R")
source("GAPIT.Numericalization.Modified.R")
setwd(home.dir)

#Read in the appropriate files
tabular.summary <- read.table(paste(home.dir,"\\(31) Tabular Summary Info for 2015 analysis\\LOD scores\\Tab_Sum_Carot_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions_LODscores.txt", sep = ""), head = TRUE)
absolute.final.data.set.FPKM <- read.table(paste(home.dir,FPKM.file.dir,"FPKM.table.by.gene.annotation.complete.founder.matrix.txt", sep = ""), head = TRUE)
										
		
###generate list of QTL numbers and common support intervals for GMMS4 script
all.JL.QTL.SI <- read.table("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Common_SI_array_for_tri_auto_June_2015.txt", head=TRUE)
			
			
### automating generation of triangulation files from all support intervals
###	using the loops from "extracting_imputed_genotype_matrix_automate.r" script		
for (q in 16:39){
  common.SI <- all.JL.QTL.SI[q,1]
  gene <- as.character(all.JL.QTL.SI[q,5])

  dir.create(paste("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_R\\QTL_",common.SI,"_imputed.ordered.tri.files.for.",gene, sep=""))
  output.folder <- (paste("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_R\\QTL_",common.SI,"_imputed.ordered.tri.files.for.",gene, "\\", sep=""))


  ### regenerate loop for trait list
  trait.list <- NULL
  for (m in 6:13){    #columns with significant traits
    select.trait <- as.character(all.JL.QTL.SI[q,m])
    trait.list <- c(trait.list, select.trait)  
  } #end m traits
  trait.list <- trait.list[!is.na(trait.list)]			

  ### regenerate loop to obtain all significant markers for these traits
  peak.marker.list <- NULL
  for (n in 14:21){    #columns with significant traits
    select.marker <- as.character(all.JL.QTL.SI[q,n])
    peak.marker.list <- c(peak.marker.list, select.marker)  
    
  } #end n peak markers in k interval
  peak.marker.list <- peak.marker.list[!is.na(peak.marker.list)]	
		
#print(paste("Following_parameters_being_used:QTL_",QTL,";gene_",gene,";chr_",chr,";trait_",trait.list,"_with_peak_marker_",peak.marker.list,sep = ""))	
 ### note, peak.marker.list must pair with traits in trait.list order
		         	
#generate matrix of traits and markers for each QTL
list.by.QTL <- cbind(trait.list, peak.marker.list)

for (p in 1:nrow(list.by.QTL)){
trait <- as.character(list.by.QTL[p,1])
QTL <- as.character(list.by.QTL[p,2])


#Obtain chromosome, start and stop support interval positions, from the "tabular summary"
row.of.interest <- which((tabular.summary[,1] == trait)  & (tabular.summary[,6] == QTL))
chr <- tabular.summary[row.of.interest, 2]
QTL.start.bp <- tabular.summary[row.of.interest, 4]
QTL.stop.bp <- tabular.summary[row.of.interest, 5]

print(paste("Following_parameters_being_used:QTL_",common.SI,";gene_",gene,";chr_",chr,";trait_",trait,"_with_peak_marker_",QTL,sep = ""))

#add an object that has all of the names of the NAM founders
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

effect.estimates.from.entire.trait <- read.table(paste(home.dir,location.of.effect.estimates,
                                      "Pop.by.Marker.Effect.Estimates.from.R.",trait,".SI01_2015.txt", sep = ""),head = TRUE)#} 
                                      ####IMPORTANT!!!!!!!!!!!  Requirement for transformed estimates only changed as of June 2015


markerID <-  substr(effect.estimates.from.entire.trait[,1], 7,10000)
popID <- substr(effect.estimates.from.entire.trait[,1], 1,5)

#transformed.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
                                #as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), ncol(effect.estimates.from.entire.trait)])

transformed.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
                                as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), 2])
                                )
colnames(transformed.effect.estimates.from.QTL) <- c("Population","Trans.Est")

transformed.effect.estimates.from.QTL <- merge(transformed.effect.estimates.from.QTL, NAM.pops, by.x = "Population", by.y = "Pop.num")


#########################################################################
#Run get me my SNPs to obtain the genotypes of these SNPs, read them in below
imputed.file.location <- paste("QTL_",common.SI,"_imputed_matrix_for_",gene,"\\", sep = "")
setwd(paste(home.dir,get.me.my.SNPs.files.dir,imputed.file.location,sep = ""))
###GWAS.SNPs.on.Chr <- GAPIT.Fragment(file.G=paste("Output_GWAS_Results_",trait,"_lt_0.05_Chr_",sep = ""),
                     ### file.Ext.G= "txt", file = chr, frag = 1, genoFormat="hapmap", file.fragment=300)


### In cases where there are no GWAS SNPs in interval and no imputed file was generated, use a mock monomorphic marker file as a placeholder for correlations
### When looking at results, will need to identify which results were derived from mock file and delete them   
   
file.name <- paste("imputed_RMIP_SNPs_",gene,"_interval_",trait,"_QTL",common.SI,"_chr",chr,".txt",sep = "")   
         
 if(file.exists(file.name)){
   real.file <- read.delim(file.name)
   setwd(paste(home.dir,get.me.my.SNPs.files.dir,imputed.file.location,sep = ""))
   GWAS.SNPs.on.Chr <- GAPIT.Fragment(file.G=paste("imputed_RMIP_SNPs_",gene,"_interval_",trait,"_QTL",common.SI,"_chr",sep = ""),                       
                      file.Ext.G= "txt", file = chr, frag = 1, genoFormat="hapmap", file.fragment=300)
   }else{
   setwd(paste(home.dir,get.me.my.SNPs.files.dir,sep = "")) 
   GWAS.SNPs.on.Chr <- GAPIT.Fragment(file.G=paste("mock.placeholder.hapmap.snp.file.chr",sep = ""),                       #######################
                      file.Ext.G= "txt", file = chr, frag = 1, genoFormat="hapmap", file.fragment=300) 
   setwd(paste(home.dir,get.me.my.SNPs.files.dir,imputed.file.location,sep = ""))  
   print(paste("-------------------- Mock_marker_file_used_for_SNP_correlations_for_QTL_",common.SI,"_gene_",gene,"_trait_",trait,"-------------------",sep = "")) 
   }   
  
  
numeric.GWAS.SNPs.data.temp <- as.data.frame(t(GWAS.SNPs.on.Chr$GD[-c(1:11)]))
colnames(numeric.GWAS.SNPs.data.temp) <- GWAS.SNPs.on.Chr$GD[,1]
final.labels.for.output.table <- paste("Chr_",GWAS.SNPs.on.Chr$GD[,3],"_",GWAS.SNPs.on.Chr$GD[,4], "_bp", sep = "")
setwd(home.dir)

############

#Obtain a vector of all of the GRZM ids of the genes within the interval
condition.1 <- absolute.final.data.set.FPKM[,3] == paste("chr", chr,sep = "")
condition.2 <- absolute.final.data.set.FPKM[,4] >= QTL.start.bp
condition.3 <- absolute.final.data.set.FPKM[,5] <= QTL.stop.bp

gene.candidate.list.from.SI <- absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3), c(1,2,4,5)]     #used later for ordering effect est data                 

gene.candidate.list.from.SI <- gene.candidate.list.from.SI[order(gene.candidate.list.from.SI[,3]),]
list.of.gene.GRZMs.in.support.interval <- as.character(gene.candidate.list.from.SI[,1])                                                      
list.of.gene.names <- as.character(gene.candidate.list.from.SI[,2])                                                                           
   

#Correlate each SNP with the JL QTL effect estimates
JL.Correlation.vector <- NULL
for(i in 1:ncol(numeric.GWAS.SNPs.data.temp)){
  #Merge the JL effect estiamtes and the ith SNP into one file
  the.SNP.of.interest <- data.frame(cbind(rownames(numeric.GWAS.SNPs.data.temp), numeric.GWAS.SNPs.data.temp[,i]))
  JL.and.SNP.alleles <- merge(transformed.effect.estimates.from.QTL, the.SNP.of.interest, by.x = "Pop.Founders", by.y = "X1")
  
  #Calcualte the Spearman rank correlation coefficient
  Correl.Spearman <- cor(JL.and.SNP.alleles[,3], as.numeric(as.character(JL.and.SNP.alleles[,4])), use = "pairwise.complete.obs", method = "spearman")
  
  #Append it to JL.Correlation.vector
  JL.Correlation.vector <- c(JL.Correlation.vector, Correl.Spearman)

}#end for(i in 1:ncol(numeric.GWAS.SNPs.data.temp))
final.output.table <- JL.Correlation.vector

#Append the name of the QTL onto final.output.table
final.output.table <- c(paste("Chr_",chr,"_QTL_Centered_at_",QTL,sep = ""),final.output.table)

#########################



for(i in list.of.gene.GRZMs.in.support.interval){
  #Correlate FPKM from each time point with JL QTL effect estimates. Use Generating_uniform_FPKM_matrix.r in this step
  the.FPKM.values <- get.me.my.FPKM.values(absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, 
                                         gene.ID = i, gene.name = list.of.gene.names[which(list.of.gene.GRZMs.in.support.interval == i)], 
                                         print.out.results = FALSE, output.dir = output.dir, home.dir = home.dir)

  #Correlate FPKM from each time point with each GWAS SNP. Highlight the SNPs that are within the JL support interval.

  for(j in 1:ncol(the.FPKM.values)){
    the.FPKM.values.of.interest <- data.frame(cbind(rownames(the.FPKM.values), the.FPKM.values[,j]))
    #Loop through the columns of the columns of numeric.GWAS.SNPs.data.temp
    correl.vector.temp <- NULL
    for(k in 1:ncol(numeric.GWAS.SNPs.data.temp)){
        #Merge the JL effect estiamtes and the ith SNP into one file
        the.SNP.of.interest <- data.frame(cbind(rownames(numeric.GWAS.SNPs.data.temp), numeric.GWAS.SNPs.data.temp[,k]))
        SNP.alleles.and.FPKM <- merge(the.SNP.of.interest, the.FPKM.values.of.interest, by.x = "X1", by.y = "X1")
        
        #Calcualte the Spearman rank correlation coefficient
        Correl.Spearman <- cor(as.numeric(as.character(SNP.alleles.and.FPKM[,2])), as.numeric(as.character(SNP.alleles.and.FPKM[,3])), 
                              use = "pairwise.complete.obs", method = "spearman")
        correl.vector.temp <- c(correl.vector.temp, Correl.Spearman)
    }#end for(k in 1:ncol(numeric.GWAS.SNPs.data.temp))
    correl.vector.temp <- c(paste(i,"_",colnames(the.FPKM.values)[j],sep = ""), correl.vector.temp)
    final.output.table <- rbind(final.output.table, correl.vector.temp)
  }#end for(j in 1:ncol(the.FPKM.values))

}#End for for(i in list.of.gene.GRZMs.in.support.interval)

#Set the column names of final output table equal to final.labels.for.output.table
colnames(final.output.table) <- c("Gene_in_SI_x_FPKM_Sample", paste("Chr_",GWAS.SNPs.on.Chr$GD[,3],"_Pos_",GWAS.SNPs.on.Chr$GD[,4],
                                  sep = ""))

#Write the results to an output table
write.table(final.output.table, paste(output.folder, "Correl_Between_FPKM_JL_GWAS_",trait, "_",QTL,"_IMPUTED_2015.txt", sep = ""),
                                  quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)


#Correlate the FPKM values with the allelic effect estimates. The output file will have the genes in the rows, and the time points
# in the columns

count <- 0
for(i in list.of.gene.GRZMs.in.support.interval){
  #Correlate FPKM from each time point with JL QTL effect estimates. Use Generating_uniform_FPKM_matrix.r in this step
  the.FPKM.values <- get.me.my.FPKM.values(absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, 
                                         gene.ID = i, gene.name = list.of.gene.names[which(list.of.gene.GRZMs.in.support.interval == i)], 
                                         print.out.results = FALSE, output.dir = output.dir, home.dir = home.dir)
  
  temp.result.vector <- i                                       
  for(j in 1:ncol(the.FPKM.values)){
    #get the jth column of the.FPKM.values
    this.FPKM <- data.frame(cbind(rownames(the.FPKM.values), the.FPKM.values[,j]))
    
    #merge it with the allelic effect estimates
    JL.effects.and.this.FPKM <- merge(transformed.effect.estimates.from.QTL, this.FPKM, by.x = "Pop.Founders", by.y = "X1")
    
    #calculate Spearman's rank correlation coefficient
    Correl.Spearman <- cor(as.numeric(as.character(JL.effects.and.this.FPKM[,3])), as.numeric(as.character(JL.effects.and.this.FPKM[,4])), 
                          use = "pairwise.complete.obs", method = "spearman")    
    
    #append it to a vector of results
    temp.result.vector <- c(temp.result.vector, Correl.Spearman)
    
  }#end for(j in 1:ncol(the.FPKM.values))
  
  
  #transformed.effect.estimates.from.QTL
  count <- count+1
  
  #Append results to the final.output.FPKM.QTL.effect.table 
  if(count == 1){
    final.output.FPKM.QTL.effect.table <- temp.result.vector
  }else{
    final.output.FPKM.QTL.effect.table <- rbind(final.output.FPKM.QTL.effect.table, temp.result.vector)
  }#end if(count == 1)
  
  
}#end for(i in list.of.gene.GRZMs.in.support.interval) 
colnames(final.output.FPKM.QTL.effect.table) <- c("GRZM_ID", colnames(the.FPKM.values))


#Order results according to genomic position
final.output.FPKM.QTL.effect.table <- as.data.frame(final.output.FPKM.QTL.effect.table)
order.final <- merge(gene.candidate.list.from.SI, final.output.FPKM.QTL.effect.table, by.x = "gene_locus", by.y = "GRZM_ID")
final.output.FPKM.QTL.effect.table.ordered <- order.final[order(order.final[,3]),]

#Export final.output.FPKM.QTL.effect.table
write.table(final.output.FPKM.QTL.effect.table.ordered, paste(output.folder, "Correl_Between_FPKM_and_",trait, "_",QTL,"_Eff_Ests_ORDERED_2015.txt", sep = ""),
                                  quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)


 } # end p trait/ marker combo from list.by.QTL
} #end q QTL intervals

