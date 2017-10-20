rm(list = ls())

###Required Files:
### (1) FPKM matrix with log2 transformed data, organized by founder, named FPKM.table.by.gene.annotation.complete.founder.matrix.txt                        
### (2) GAPIT source files                                   
### (3) Tabular summary file                                
### (4) TRANSFORMED effect estimate files, from script (2-2)   


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
    setwd(output.dir)
    write.table(FPKM.matrix, paste("FPKM.matrix.for.",gene.name, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
  }
  return(FPKM.matrix)
}# end get.me.my.FPKM.values

#########################################################################
correlate.it.dude <- function(transformed.effect.estimates.from.QTL =NULL, tabular.summary = tabular.summary, absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, trait = trait, QTL=QTL){
  #Obtain chromosome, start and stop support interval positions, from the "tabular summary"
  row.of.interest <- which((tabular.summary[,1] == trait)  & (tabular.summary[,6] == QTL))
  chr <- tabular.summary[row.of.interest, 2]
  QTL.start.bp <- tabular.summary[row.of.interest, 4]
  QTL.stop.bp <- tabular.summary[row.of.interest, 5]
   
  #Obtain a vector of all of the GRZM ids of the genes within the interval
  condition.1 <- absolute.final.data.set.FPKM[,3] == paste("chr", chr,sep = "")
  condition.2 <- absolute.final.data.set.FPKM[,4] >= QTL.start.bp
  condition.3 <- absolute.final.data.set.FPKM[,5] <= QTL.stop.bp
  
  list.of.gene.GRZMs.in.support.interval <- as.character(absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3),1])
  list.of.gene.names <- as.character(absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3),2])
  list.of.chr <- as.character(absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3),3])
  list.of.bp.start <- absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3),4]
  list.of.bp.stop <- absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3),5]
   
  print("-------------------------Correlating each FPKM values with allelic effect estimates--------------------------------------------------------------")
  #Correlate the FPKM values with the allelic effect estimates. The output file will have the genes in the rows, and the time points in the columns
  
  count <- 0
  for(i in list.of.gene.GRZMs.in.support.interval){
    #Correlate FPKM from each time point with JL QTL effect estimates. Use Generating_uniform_FPKM_matrix.r in this step
    the.FPKM.values <- get.me.my.FPKM.values(absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, 
                                             gene.ID = i, gene.name = list.of.gene.names[which(list.of.gene.GRZMs.in.support.interval == i)], 
                                             print.out.results = FALSE, output.dir = output.dir, home.dir = home.dir)
    temp.result.vector.Spearman <- i   
    temp.result.vector.Pearson <- i                               
    for(j in 1:ncol(the.FPKM.values)){
      #get the jth column of the.FPKM.values
      this.FPKM <- data.frame(cbind(rownames(the.FPKM.values), the.FPKM.values[,j]))
      
      #merge it with the allelic effect estimates
      JL.effects.and.this.FPKM <- merge(transformed.effect.estimates.from.QTL, this.FPKM, by.x = "Pop.Founders", by.y = "X1")
      
      #calculate Spearman and Pearson - NOTE: USING PEARSON IN FINAL OUTPUT
      Correl.Spearman <- suppressWarnings(cor(as.numeric(as.character(JL.effects.and.this.FPKM[,3])), as.numeric(as.character(JL.effects.and.this.FPKM[,4])), 
                                              use = "pairwise.complete.obs", method = "spearman"))
      tryCatch({
        Correl.Pearson <- cor(as.numeric(as.character(JL.effects.and.this.FPKM[,3])), as.numeric(as.character(JL.effects.and.this.FPKM[,4])), 
            use = "pairwise.complete.obs", method = "pearson") 
      }, warning = function(w) {
        Correl.Pearson <- NA
        these.effects = as.vector(JL.effects.and.this.FPKM[,3])
        these.FPKMs = as.vector(JL.effects.and.this.FPKM[,4])
        print(paste(QTL,trait,i,", time point in col",j,sep=" "))
        print(paste(c("Produced NA for FPKM/JL SNP allelic effect correlation because stdev was 0. FPKM vector ",these.FPKMs," Effect vector ",these.effects),collapse=','))
        }
      )
      
      #append to results vector
      temp.result.vector.Spearman <- c(temp.result.vector.Spearman, Correl.Spearman)
      temp.result.vector.Pearson <- c(temp.result.vector.Pearson, Correl.Pearson)
     
    }#end for(j in 1:ncol(the.FPKM.values))
    #transformed.effect.estimates.from.QTL
    count <- count+1
    
    #Append results to the final.output.FPKM.QTL.effect.table 
    if(count == 1){
      final.output.FPKM.QTL.effect.table <- temp.result.vector.Pearson
    }else{
      final.output.FPKM.QTL.effect.table <- rbind(final.output.FPKM.QTL.effect.table, temp.result.vector.Pearson)
    }#end if(count == 1)
  }#end for(i in list.of.gene.GRZMs.in.support.interval) 
  
  #Put the column names on final.output.FPKM.QTL.effect.table
  colnames(final.output.FPKM.QTL.effect.table) <- c("GRZM_ID", colnames(the.FPKM.values))

  #Merge the the information about chromosome and bp position to "final.output.FPKM.QTL.effect.table"
  final.output.FPKM.QTL.effect.table <- as.data.frame(final.output.FPKM.QTL.effect.table)
  data.frame.with.chr.pos <- as.data.frame(cbind(list.of.gene.GRZMs.in.support.interval, list.of.gene.names, list.of.chr, list.of.bp.start, list.of.bp.stop, rep(QTL, length(list.of.bp.start))) ) 
  colnames(data.frame.with.chr.pos) <- c("GRZM_ID", "Gene_Name", "Chr", "Start_bp", "Stop_bp", "JL_QTL")
  
  final.output.FPKM.QTL.effect.table <- merge(final.output.FPKM.QTL.effect.table, data.frame.with.chr.pos, by.x = "GRZM_ID", by.y = "GRZM_ID")

  #Sort this table in genomic order
  final.output.FPKM.QTL.effect.table <- final.output.FPKM.QTL.effect.table[order(final.output.FPKM.QTL.effect.table[,12]),] 
  
  return(final.output.FPKM.QTL.effect.table)
} #end correlate.it.dude()

#####Obtain the P-value for testing H0:rho = 0
obtain.P.value <- function(correl = NA,i=NULL){
  #dfr <- 20 # = 22 - 2 (Number of observations with "good" FPKM measurements minus 2; we loose two d.f.; one for estimating the mean of each variable)
  dfr = nfam_perTimePt[i-1] #subtract 1 because first time point has i index of 2
  print(paste("dfr actually passed to function is ",dfr,sep=""))
  r2 <- correl^2
  Fstat <- r2 * dfr / (1 - r2)
  P.val <- 1 - pf(Fstat, 1, dfr) 
  return(P.val)
}#end obtain.P.value()

GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure.mod <-
  function(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH"){
    #Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
    #Output: PWIP, number.of.significant.SNPs
    #Authors: Alex Lipka and Zhiwu Zhang 
    # Last update: May 5, 2011 
    ##############################################################################################
    #Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest
    if(is.null(PWI))
    {
      PWIP=NULL
      number.of.significant.SNPs = 0
    }
    
    if(!is.null(PWI))
    {  
      if(length(which(is.na(PWI)))>0){
        PWI[which(is.na(PWI))] = 1
      }

    for(i in 1:ncol(PWI))
    {
      #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
      #Right column: Simes corrected p-value
      res <- mt.rawp2adjp(PWI[,i], FDR.Procedure)
      
      #This command should order the p-values in the order of the SNPs in the data set
      adjp <- res$adjp[order(res$index), ]
      
      #round(adjp[1:7,],4)
      #Logical statment: 0, if Ho is not rejected; 1, if  Ho is rejected, by the Simes corrected p-value
      #  temp <- mt.reject(adjp[,2], FDR.Rate)
      
      #Lists all number of SNPs that were rejected by the BY procedure
      #temp$r
      
      #Attach the FDR adjusted p-values to AS_Results
      if(i == 1){
        PWIP <- adjp[,2]
      }else{
        PWIP <- cbind(PWIP, adjp[,2])
      }

      #Sort these data by lowest to highest FDR adjusted p-value
      #PWIP <- PWIP[order(PWIP[,4]),]
      
    } #end for(i in 1:ncol(PWI))
    
    }
    #return(list(PWIP=PWIP, number.of.significant.SNPs = number.of.significant.SNPs))
    return(PWIP)
    
  }#end GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure.mod()

get.me.my.results <- function(tabular.summary = tabular.summary, absolute.final.data.set.FPKM, trait = trait){ 
  
  the.QTL.for.the.trait <- as.vector(as.character(tabular.summary[which(as.character(tabular.summary[,1]) == trait),6]))
  print(paste(the.QTL.for.the.trait,sep = ""))
  the.correlation.results <- NULL
  for(QTL in the.QTL.for.the.trait){
    #QTL <- QTL.under.study
    #Gene <- NA
    #Gene.GRZM.ID <- NA
    print(paste("QTL = ", QTL,sep = ""))
    transformed.effect.estimates.from.QTL <- the.extractinator(tabular.summary = tabular.summary, QTL = QTL, Gene =Gene, Gene.GRZM.ID = Gene.GRZM.ID, trait = trait)
    results.summaries.QTL.under.study <- correlate.it.dude(transformed.effect.estimates.from.QTL = transformed.effect.estimates.from.QTL, QTL=QTL,
                                                           tabular.summary = tabular.summary, absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, trait = trait)
    the.correlation.results <- rbind(the.correlation.results, results.summaries.QTL.under.study)
    print(dim(the.correlation.results))
  }#end for(QTL.under.study in the.QTL.for.the.trait)
  print("---Dimension of the.correlation.results after the loop---")
  print(dim(the.correlation.results))
  #Obtain the P-values for testing H0:rho = 0 
  for(i in 2:9){
    temp.vector <- as.numeric(as.vector(the.correlation.results[,i]))
    vector.of.P.values <- unlist(lapply(temp.vector, obtain.P.value,i=i))
    if(i == 2){
      matrix.of.P.values <- vector.of.P.values
    }else{
      matrix.of.P.values <- cbind(matrix.of.P.values, vector.of.P.values)
    }
  }
  the.P.value.results <- cbind(the.correlation.results[,1], matrix.of.P.values, the.correlation.results[,10:14])
  colnames(the.P.value.results)[2:9] <- colnames(the.correlation.results)[2:9]
  colnames(the.P.value.results)[1] <-colnames(the.correlation.results)[1]
  
  #Obtain the FDR-adjusted P-values
  the.FDR.Adjusted.P.values <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure.mod(PWI = matrix.of.P.values)
  
  print(paste("----------------And now, we are calculating the FDR-adjusted P-values for ", trait, "------------------------", sep = ""))
  
  the.FDR.adj.P.value.results <- cbind(the.correlation.results[,1], the.FDR.Adjusted.P.values, the.correlation.results[,10:14])
  colnames(the.FDR.adj.P.value.results)[2:9] <- colnames(the.correlation.results)[2:9]
  colnames(the.FDR.adj.P.value.results)[1] <-colnames(the.correlation.results)[1]
  print("-------Results are imminent-----------")
  bonferroni.adjusted.alpha <- 0.05/6
  return(list(the.correlation.results = the.correlation.results, the.P.value.results = the.P.value.results, 
              the.FDR.adj.P.value.results = the.FDR.adj.P.value.results, bonferroni.adjusted.alpha=bonferroni.adjusted.alpha))
  
}#end get.me.my.results() 


###!NOTE: The "get me my SNPs" output file HAS to be labeled as follows, e.g. "Output_GWAS_Results_AT_lt_0.05_Chr_5.txt"
the.extractinator <- function(tabular.summary = tabular.summary, QTL = QTL, Gene = NA, Gene.GRZM.ID = NA, trait = trait){
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
  #7/21/2015 CHD changed Il14H to IL14H to match nomenclature above
  founder.names <- as.data.frame(c("B97", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", 
                                   "CML69", "HP301", "IL14H", "KI11", "KI3", "KY21", "M162W", "M37W", "MO18W", 
                                   "MS71", "NC350", "NC358", "OH43", "OH7B", "P39", "TX303", "TZI8"))
  NAM.pops <- cbind(pop.seq, founder.names)
  colnames(NAM.pops) <- c("Pop.num", "Pop.Founders")
  
  #Obtain the populationxmarker effect estimates for this QTL
  if(trait.set == "Carot"){
  effect.estimates.from.entire.trait <- read.table(paste(location.of.effect.estimates,
                                                         "Pop.by.Marker.Effect.Estimates.from.R.",trait,".SI01_2015.txt", sep = ""),head = TRUE)  
  }  
  if(trait.set == "tocos"){
    if(trait=="dT3"){
      effect.estimates.from.entire.trait <- read.table(paste(dT3.new.path,
                                                             "Pop.by.Marker.Effect.Estimates.from.R.dT3_Redone.SI01.txt", sep = ""),head = TRUE)#} 
    }else{
      effect.estimates.from.entire.trait <- read.table(paste(location.of.effect.estimates,
                                                             "Pop.by.Marker.Effect.Estimates.from.R.",trait,".SI01.txt", sep = ""),head = TRUE)#} 
      ####IMPORTANT!!!!!!!!!!!  Requirement for transformed estimates only changed as of June 2015
    }
  }
  
  markerID <-  substr(effect.estimates.from.entire.trait[,1], 7,10000)
  popID <- substr(effect.estimates.from.entire.trait[,1], 1,5)
  
  ##The below chunk is from before lambda sign correction was implemented.
  # #transformed.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
  # #                                          as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), ncol(effect.estimates.from.entire.trait)])
  # #)
  # 
  
  transformed.effect.estimates.from.QTL = NULL #to make sure cleaned out from previous iteration
  
  transformed.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
                                as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), 2])
                                )
  
  #This part was brought over from 5-2 Pearson lambda script, 9/20/2017
  lambda = lambda.values[which(lambda.values[,2] == trait),3]
  scalar = if(lambda < 0){-1}else{1}  
  
  print("testing effect signs BEFORE multiplying by -1 if lambda negative.")
  print(transformed.effect.estimates.from.QTL[1:10,2])
  transformed.effect.estimates.from.QTL[,2] = transformed.effect.estimates.from.QTL[,2]*scalar
  print("testing effect signs AFTER multiplying by -1 if lambda negative.")
  print(transformed.effect.estimates.from.QTL[1:10,2])
  
  colnames(transformed.effect.estimates.from.QTL) <- c("Population","Trans.Est")
  
  transformed.effect.estimates.from.QTL <- merge(transformed.effect.estimates.from.QTL, NAM.pops, by.x = "Population", by.y = "Pop.num")

  return(transformed.effect.estimates.from.QTL)
}#end the.extractinator 


##########################################################################
##########################################################################
##########################################################################


library(multtest)
trait.set = "tocos" #Options are "Carot" or "tocos"
FDR_corr = "by_QTL" #Options are "genomewide" or "by_QTL"
#ran genomewide on 30 July 2016 to get corrected raw P-values, then by_QTL on 31 July to obtain corrected FDR by Q(correct nFAM)
if(trait.set == "Carot"){
  setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
  home.dir <- getwd()
  location.of.modified.GAPIT.files <- "\\(21) GAPIT source files\\"
  location.of.effect.estimates <- "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\"
  location.of.GWAS.results <- "\\(10)GWAS Analysis\\RUV GWAS 25fam_alldata_alpha01_2015_FINAL\\"
  output.dir <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Significance Threshold\\FDR_Test_Ancillary_Output_Files\\"
  get.me.my.SNPs.files.dir <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Results_from_GMMS4_imputed_matrix_T3_2015\\" 
  FPKM.file.dir <- "\\(15)Correlated Expression\\"
  correlation.output.dir <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015\\Significance Threshold\\FDR_Corrected_Corr_Results\\"
  tabSummary.path = "\\(31) Tabular Summary Info for 2015 analysis\\LOD scores\\"
}

if(trait.set == "tocos"){
  setwd("C:/Users/chd45/Documents/Data_From_Mac/Documents/Toco_NAM/")
  home.dir <- getwd()
  tabSummary.path = paste(home.dir,"/Summary_Tables.Figures/",sep='')
  location.of.modified.GAPIT.files <- paste(home.dir,"/Expression/inputsFor5-2/",sep='')
  location.of.effect.estimates <- paste(home.dir,"/JL/Allelic_Effect_Estimates.no.MultiColl/",sep='')
  dT3.new.path = paste(home.dir,"/Methods/dT3_removeExtremeVal_test/new.trans_new.perm_FINAL/",sep='')
  location.of.GWAS.results <- paste(home.dir,"/GWAS/Summaries_by_Trait/",sep='')
  #dir.for.5.2.ordered <- paste(home.dir,"/Expression/orderedFrom5.2/",sep='')
  dir.for.5.1.matrices <- paste(home.dir,"/Expression/ImputedMatricesfrom5.1/",sep='')
  FPKM.file.dir <- "C:/Users/chd45/Documents/Data_From_Mac/Documents/Toco_NAM/for_FPKMxEffect_comparison/"
  FPKM.output.dir = paste(home.dir,"/Expression/FPKMmatricesFrom5.3/",sep='')
  correlation.output.dir <- paste(home.dir,"/Expression/corrsFrom5.3/",sep='')
  fdr.adjusted.for.each.SI.dir <- correlation.output.dir
  lambda.dir = "C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/"
}

tabular.summary <- read.table(paste(tabSummary.path,"Tab_Sum_Final_dT3_Redone_with_Common_SI_Info_left_bound.txt",sep=''), head = TRUE,stringsAsFactors=FALSE)
trait.list <- as.character(unique(tabular.summary[,1]))
lambda.values <-read.table(paste(lambda.dir,"BLUEs_lambda_values_used_all_dT3.only.Redone.txt", sep = ""), head = TRUE)
#trait.list = "dT3" #for re-do

#Source in the modified GAPIT files
setwd(location.of.modified.GAPIT.files)
source("GAPIT.Fragment.Modified.R")
source("GAPIT.HapMap.Modified.R")
source("GAPIT.Numericalization.Modified.R")
setwd(home.dir)

#Read in the appropriate files
absolute.final.data.set.FPKM <- read.table(paste(FPKM.file.dir,"FPKM.table.by.gene.ann.complete.matrix.FPKMthreshold.1_filter.by.kernel_across_all.samples.log2trans.txt", sep = ""), head = TRUE)
nfam_perTimePt = c(21,21,19,20,20,19,20,20) #number of families for each time point, in same order iterated through in get.me.my.results (12,16,20,24,30,36,101,202)
#nfam_perTimePt = c(19,19,17,18,18,17,18,18) #old INCORRECT #s; B73 and Mo17 were subtracted twice

#Get the results we need for both traits
if(FDR_corr == "genomewide"){
  for(trait in trait.list){
    results.test <- get.me.my.results(tabular.summary = tabular.summary, absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, trait = trait)
    write.table(results.test$the.correlation.results, paste(correlation.output.dir,"Correlations.for.",trait,"_Correct.nFam_Pearson_lambda.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
    write.table(results.test$the.P.value.results, paste(correlation.output.dir,"Raw.P.values.for.",trait,"_Correct.nFam_Pearson_lambda.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
    write.table(results.test$the.FDR.adj.P.value.results, paste(correlation.output.dir,"FDR.Adjusted.P.values.for.",trait,"_Correct.nFam_Pearson_lambda.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
}#end for(trait in unique(tabular.summary[,1]))
}#end if genomewide

if(FDR_corr == "by_QTL"){
  for(trait in trait.list){
  setwd(correlation.output.dir)
  the.raw.P.values <- read.table(paste("Raw.P.values.for.", trait, "_Correct.nFam_Pearson_lambda.txt", sep = ""),head = TRUE)
  
  #For each unique QTL
  count <- 0
  
  for(QTL in unique(the.raw.P.values[,ncol(the.raw.P.values)])){
    #Split up the data so that you only get the results for the QTL under study
    the.raw.P.values.for.this.QTL <- the.raw.P.values[which(the.raw.P.values[,ncol(the.raw.P.values)] == QTL),]
    
    #Run the FDR correction
    the.FDR.Adjusted.P.values <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure.mod(PWI = the.raw.P.values.for.this.QTL[,2:9])
    
    print(paste("----------------And now, we are calculating the FDR-adjusted P-values BY QTL for ", trait, "------------------------", sep = ""))
    
    the.FDR.adj.P.value.results <- cbind(the.raw.P.values.for.this.QTL[,1], the.FDR.Adjusted.P.values, 
                                         the.raw.P.values.for.this.QTL[,10:14],nrow(the.raw.P.values.for.this.QTL))
    colnames(the.FDR.adj.P.value.results)[2:9] <- colnames(the.raw.P.values.for.this.QTL)[2:9]
    colnames(the.FDR.adj.P.value.results)[1] <-colnames(the.raw.P.values.for.this.QTL)[1]
    colnames(the.FDR.adj.P.value.results)[15] <- "nGenes_this.cSI_thisTrait"
    
    #Append the support interval-wise FDR corrected P-values
    if(count == 0){
      FDR.adjusted.P.values.for.all.SIs <- the.FDR.adj.P.value.results 
    }else{
      FDR.adjusted.P.values.for.all.SIs <- rbind(FDR.adjusted.P.values.for.all.SIs,the.FDR.adj.P.value.results) 
    }
    count <- count + 1
  }#End for each unique QTL
  
  #Export the FDR-adjusted P-values, adjusted separately for each support interval
  setwd(fdr.adjusted.for.each.SI.dir)
  write.table(FDR.adjusted.P.values.for.all.SIs, paste("FDR.Adjusted.P.values.for.",trait,"_by_Q_Correct.nFam_Pearson_lambda.txt", sep = ""),
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  }#end if output by QTL
}#end for trait in trait.list
