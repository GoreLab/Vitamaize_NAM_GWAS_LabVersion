rm(list = ls())

###Required Files:
### (1) FPKM matrix with log2 transformed data, organized by founder, named FPKM.table.by.gene.annotation.complete.founder.matrix.txt
### (2) GAPIT source files
### (3) Tabular summary file
### (4) TRANSFORMED effect estimate files, from script (2-2)
### (5) GWAS Results Summary files, from script (3-1)
### (6) QTL array file called in script (5-1)
### (7) Imputed and extracted genotype matrices from script (5-1)


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

##################################################################################################################################################################

correlate.GWAS.SNPs.with.JL.Effect.Estims = function(numeric.GWAS.SNPs.data.temp = NA, transformed.effect.estimates.from.QTL = NA){
  #Correlate each SNP with the JL QTL effect estimates
  print(paste("Now correlating GWAS SNPs with JL QTL effect estimates for trait",trait," and common SI",q,sep=''))
  JL.Correlation.vector <- NULL
  for(i in 1:ncol(numeric.GWAS.SNPs.data.temp)){
    #Merge the JL effect estiamtes and the ith SNP into one file
    the.SNP.of.interest <- data.frame(cbind(rownames(numeric.GWAS.SNPs.data.temp), numeric.GWAS.SNPs.data.temp[,i]))
    JL.and.SNP.alleles <- merge(transformed.effect.estimates.from.QTL, the.SNP.of.interest, by.x = "Pop.Founders", by.y = "X1")
  
    #Calculate Spearman rank and Pearson coefficients - NOTE: ONLY APPENDING PEARSON TO FINAL OUTPUT
    Correl.Spearman <- suppressWarnings(cor(JL.and.SNP.alleles[,3], as.numeric(as.character(JL.and.SNP.alleles[,4])), use = "pairwise.complete.obs", method = "spearman"))
    tryCatch({
      Correl.Pearson <- cor(JL.and.SNP.alleles[,3], as.numeric(as.character(JL.and.SNP.alleles[,4])), use = "pairwise.complete.obs", method = "pearson")
    }, warning = function(w) {
      Correl.Pearson = NA
      these.effects = as.vector(JL.and.SNP.alleles[,3])
      these.alleles = as.vector(JL.and.SNP.alleles[,4])
      print(paste("For common SI ",common.SI," trait ",trait, "SNP ",i,sep=""))
      print(paste(c("Produced NA for JL effect estimate/GWAS SNP correlation because stdev was 0. JL effect estimate vector ",these.effects,". GWAS SNP score vector ",these.alleles),collapse=','))
    }
    )
    
    JL.Correlation.vector <- c(JL.Correlation.vector, Correl.Pearson)  
  }#end for(i in 1:ncol(numeric.GWAS.SNPs.data.temp))
  final.output.table <- c(paste("Chr_",chr,"_QTL_Centered_at_",QTL,sep = ""),JL.Correlation.vector)
  return(final.output.table)
}#end correlate.GWAS.SNPs.with.JL.Effect.Estims

##################################################################################################################################################################

correlate.FPKMs.with.GWAS.SNPs.and.JL.Effect.Estims = function(numeric.GWAS.SNPs.data.temp = NA,transformed.effect.estimates.from.QTL = NA,final.output.table=NULL){
  print(paste("Now correlating FPKMs with GWAS SNPs for trait ",trait," and common SI ",q,sep=''))
  
  for(i in list.of.gene.GRZMs.in.support.interval){
    #Extract FPKM for all time points with JL QTL effect estimates. Use Generating_uniform_FPKM_matrix.r in this step
    the.FPKM.values <- get.me.my.FPKM.values(absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, 
                                             gene.ID = i, gene.name = list.of.gene.names[which(list.of.gene.GRZMs.in.support.interval == i)], 
                                             print.out.results = FALSE, output.dir = output.dir, home.dir = home.dir)
    
    #Correlate FPKM from each time point with each GWAS SNP. Highlight the SNPs that are within the JL support interval.
    for(j in 1:ncol(the.FPKM.values)){
      the.FPKM.values.of.interest <- data.frame(cbind(rownames(the.FPKM.values), the.FPKM.values[,j]))
      
      #Loop through columns of numeric.GWAS.SNPs.data.temp
      correl.vector.temp <- NULL
      for(k in 1:ncol(numeric.GWAS.SNPs.data.temp)){
        #Merge JL effect estimates and ith SNP into one file
        the.SNP.of.interest <- data.frame(cbind(rownames(numeric.GWAS.SNPs.data.temp), numeric.GWAS.SNPs.data.temp[,k]))
        SNP.alleles.and.FPKM <- merge(the.SNP.of.interest, the.FPKM.values.of.interest, by.x = "X1", by.y = "X1")
      
        #Calcualte Spearman rank correlation coefficient and Pearson coefficient - NOTE: USING PEARSON IN FINAL OUTPUT
        Correl.Spearman <-suppressWarnings(cor(as.numeric(as.character(SNP.alleles.and.FPKM[,2])), as.numeric(as.character(SNP.alleles.and.FPKM[,3])), 
                                           use = "pairwise.complete.obs", method = "spearman"))
        tryCatch({
          Correl.Pearson <- cor(as.numeric(as.character(SNP.alleles.and.FPKM[,2])), as.numeric(as.character(SNP.alleles.and.FPKM[,3])), 
            use = "pairwise.complete.obs", method = "pearson")
        }, warning = function(w) {
          Correl.Pearson = NA
          these.alleles = as.vector(SNP.alleles.and.FPKM[,2])
          these.FPKMs = as.vector(SNP.alleles.and.FPKM[,3])
          print(paste("For common SI ",common.SI," trait ",trait, " gene ", i ," at time point ",j, "and SNP ", k, sep=""))
          print(paste(c("Produced NA for FPKM/GWAS SNP correlation because stdev was 0. FPKM vector ",these.FPKMs,". Allele vector ",these.alleles),collapse=','))
        }
        )
      
        correl.vector.temp <- c(correl.vector.temp, Correl.Pearson)
      }#end for(k in 1:ncol(numeric.GWAS.SNPs.data.temp))
      correl.vector.temp <- c(paste(i,"_",colnames(the.FPKM.values)[j],sep = ""), correl.vector.temp)
      final.output.table <- rbind(final.output.table, correl.vector.temp)
    }#end for(j in 1:ncol(the.FPKM.values))
  }#End for(i in list.of.gene.GRZMs.in.support.interval)
  
  colnames(final.output.table) <- c("Gene_in_SI_x_FPKM_Sample", paste("Chr_",numeric.GWAS.SNPs.data.temp.with.Head[,3],"_Pos_",numeric.GWAS.SNPs.data.temp.with.Head[,4],sep = ""))
  write.table(final.output.table, paste(output.folder, "Pearson_Between_FPKM_JL_GWAS_",trait, "_",QTL,"_IMPUTED_2015.txt", sep = ""),
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  #finished calculating correlation between GWAS SNPs and JL effect estims
  
  ##################################################################################################################################################################
  #Correlate FPKM values with allelic effect estimates. output file will have genes in rows, time points in columns
  print(paste("Now correlating FPKMs with JL effect estimates for trait",trait," and common SI",q,sep=''))
  count <- 0
  for(i in list.of.gene.GRZMs.in.support.interval){
    #Correlate FPKM from each time point with JL QTL effect estimates. Use Generating_uniform_FPKM_matrix.r in this step
    the.FPKM.values <- get.me.my.FPKM.values(absolute.final.data.set.FPKM = absolute.final.data.set.FPKM, 
                                           gene.ID = i, gene.name = list.of.gene.names[which(list.of.gene.GRZMs.in.support.interval == i)], 
                                           print.out.results = FALSE, output.dir = output.dir, home.dir = home.dir)
    temp.result.vector.Spearman <- i   
    temp.result.vector.Pearson <- i                                      
    for(j in 1:ncol(the.FPKM.values)){
      this.FPKM <- data.frame(cbind(rownames(the.FPKM.values), the.FPKM.values[,j]))
      JL.effects.and.this.FPKM <- merge(transformed.effect.estimates.from.QTL, this.FPKM, by.x = "Pop.Founders", by.y = "X1")
    
      #calculate Spearman and Pearson - NOTE: USING PEARSON IN FINAL OUTPUT
      Correl.Spearman <- suppressWarnings(cor(as.numeric(as.character(JL.effects.and.this.FPKM[,3])), as.numeric(as.character(JL.effects.and.this.FPKM[,4])), 
                                            use = "pairwise.complete.obs", method = "spearman"))
      tryCatch({
        Correl.Pearson <- cor(as.numeric(as.character(JL.effects.and.this.FPKM[,3])), as.numeric(as.character(JL.effects.and.this.FPKM[,4])), 
          use = "pairwise.complete.obs", method = "pearson") 
      }, warning = function(w) {
        Correl.Pearson = NA
        these.effects = as.vector(JL.effects.and.this.FPKM[,3])
        these.FPKMs = as.vector(JL.effects.and.this.FPKM[,4])
        print(paste(common.SI,trait,i,sep=" "))
        print(paste(c("Produced NA for FPKM/JL SNP allelic effect correlation because stdev was 0. FPKM vector ",these.FPKMs," Effect vector ",these.effects),collapse=','))
      }
      )
      #append to results vector
      temp.result.vector.Spearman <- c(temp.result.vector.Spearman, Correl.Spearman)
      temp.result.vector.Pearson <- c(temp.result.vector.Pearson, Correl.Pearson)
    
    }#end for(j in 1:ncol(the.FPKM.values))
    
    count <- count+1
    if(count == 1){
      final.output.FPKM.QTL.effect.table <- temp.result.vector.Pearson
    }else{
      final.output.FPKM.QTL.effect.table <- rbind(final.output.FPKM.QTL.effect.table, temp.result.vector.Pearson)
    }#end if(count == 1)
  
  }#end for(i in list.of.gene.GRZMs.in.support.interval) 
  
  #Label results, Sort according to genomic position, Export
  colnames(final.output.FPKM.QTL.effect.table) <- c("GRZM_ID", colnames(the.FPKM.values))
  final.output.FPKM.QTL.effect.table <- as.data.frame(final.output.FPKM.QTL.effect.table)
  order.final <- merge(gene.candidate.list.from.SI, final.output.FPKM.QTL.effect.table, by.x = "gene_locus", by.y = "GRZM_ID")
  final.output.FPKM.QTL.effect.table.ordered <- order.final[order(order.final[,3]),]
  write.table(final.output.FPKM.QTL.effect.table.ordered, paste(output.folder, "Pearson_Between_FPKM_and_",trait, "_",QTL,"_Eff_Ests_ORDERED_2015.txt", sep = ""),
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  #finished calculating correlation between GWAS SNPs and FPKM expression values
  
} #end correlate.FPKMs.with.GWAS.SNPs.and.JL.Effect.Estims composite function

###################################################End of functions section; start of main code###################################################################
##################################################################################################################################################################
##################################################################################################################################################################

trait.set = "tocos" #Options are "Carot" or "tocos"

#Set the working directory
if(trait.set == "Carot"){
  setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
  home.dir <- getwd()
  tabSummary.path <- "\\(31) Tabular Summary Info for 2015 analysis\\" 
  location.of.modified.GAPIT.files <- "\\(21) GAPIT source files\\"
  location.of.effect.estimates <- "\\(9)JL Analysis\\Permutations\\Data_for_alpha01_new_TASSEL3\\Effect_estimates_TASSEL3_alpha01_2015\\"
  location.of.GWAS.results <- "\\(10)GWAS Analysis\\RUV GWAS 25fam_alldata_alpha01_2015_FINAL\\"
  dir.for.5.2.ordered <- "\\(15)Correlated Expression\\Results_from_R_FPKM1_logtyes\\"
  dir.for.5.1.matrices <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015_union\\Results_from_GMMS31_imputed_matrix_T3_2015\\"
  output.dir <- "\\(15)Correlated Expression\\Results_from_R_FPKM1_logtyes\\"
  get.me.my.SNPs.files.dir <- "\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015_union\\Results_from_GMMS31_imputed_matrix_T3_2015\\"           # this was changed
  FPKM.file.dir <- "\\(15)Correlated Expression\\"
  expression.dir <- "\\(15)Correlated Expression\\Results_from_R_FPKM1_logtyes\\"
  lambda.dir = "\\(8b)BLUE Outlier Rem and Trans Redo\\"              
  numb.common.SI = 39
  max.numb.traits.per.SI = 8
  
  tabular.summary <- read.table(paste(home.dir,tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_GWAS_FamPVE_common_SI_recsuppregions_LODscores_20150612.txt", sep = ""), head = TRUE)
  lambda.values <-read.table(paste(home.dir,lambda.dir,"Carot_BLUE_lambda_values_all.txt", sep = ""), head = TRUE)         #cbk added 7/9/15, note [,1:3] (trait.no, trait.ID, lambda)
}

if(trait.set == "tocos"){
  setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
  home.dir <- getwd()
  tabSummary.path = "/Tabular_Summaries/"
  location.of.modified.GAPIT.files <- "/Expression_Analyses/inputsFor5-2/"
  location.of.effect.estimates <- "/Allelic_Effect_Estimates.no.MultiColl/"
  location.of.GWAS.results <- "/GWAS_Analysis/GWAS_25fam_HMPonly_TASSEL3_alpha01_2015_corr/"
  dir.for.5.2.ordered <- "/Expression_Analyses/GMMS3.1Results_tocos/orderedFrom5.2/"
  dir.for.5.1.matrices <- "/Expression_Analyses/GMMS3.1Results_tocos/ImputedMatricesfrom5.1/"
  FPKM.file.dir <- location.of.modified.GAPIT.files
  expression.dir = "/Expression_Analyses/"
  lambda.dir = "C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/"
  numb.common.SI = 48
  max.numb.traits.per.SI = 9
  
  tabular.summary <- read.table(paste(home.dir,tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions_LODscores.txt", sep = ""), head = TRUE)
  lambda.values <-read.table(paste(lambda.dir,"tocos_BLUE_lambda_values_used_all.txt", sep = ""), head = TRUE) 
}

#Source in the modified GAPIT files
setwd(paste(home.dir,location.of.modified.GAPIT.files,sep = ""))
source("GAPIT.Fragment.Modified.R")
source("GAPIT.HapMap.Modified.R")
source("GAPIT.Numericalization.Modified.R")
setwd(home.dir)

#Read in the appropriate files

absolute.final.data.set.FPKM <- read.table(paste(home.dir,FPKM.file.dir,"FPKM.table.by.gene.ann.complete.matrix.FPKMthreshold.1_filter.by.kernel_across_all.samples.log2trans.txt", sep = ""), head = TRUE)

###generate list of QTL numbers and common support intervals for GMMS4 script
all.JL.QTL.SI <- read.table(paste(home.dir,expression.dir,"Common_SI_array_for_tri_auto_June_2015.txt",sep=''), head=TRUE)

sink(paste(home.dir,expression.dir,"log_Script5-2.txt",sep=''))
### automating generation of triangulation files from all support intervals
###	using the loops from "extracting_imputed_genotype_matrix_automate.r" script		
for (q in 1:numb.common.SI){
#for (q in 22:24){
  print(paste("processing Common SI number",q,sep=''))
  common.SI <- all.JL.QTL.SI[q,1]
  gene <- as.character(all.JL.QTL.SI[q,5])
  
  if(trait.set == "Carot"){
  dir.create(paste(home.dir,dir.for.5.2.ordered,"QTL_",common.SI,"_imputed.ordered.tri.files.for.",gene, sep=""),showWarnings=FALSE)
  output.folder <- (paste(home.dir,dir.for.5.2.ordered,"QTL_",common.SI,"_imputed.ordered.tri.files.for.",gene, "\\", sep=""))
  }
  if(trait.set == "tocos"){
  dir.create(paste(home.dir,dir.for.5.2.ordered,"QTL_",common.SI,"_tri.files_FPKM1_logtyes/", sep=''),showWarnings=FALSE)
  output.folder <- paste(home.dir,dir.for.5.2.ordered,"QTL_",common.SI,"_tri.files_FPKM1_logtyes/", sep='')
  }
  
  ### regenerate loop for trait list
  trait.list <- NULL
  for (m in 6:(5+max.numb.traits.per.SI)){    #columns with significant traits
    select.trait <- as.character(all.JL.QTL.SI[q,m])
    trait.list <- c(trait.list, select.trait)  
  } #end m traits
  trait.list <- trait.list[!is.na(trait.list)]			

  ### regenerate loop to obtain all significant markers for these traits
  peak.marker.list <- NULL
  for (n in (6+max.numb.traits.per.SI:(6+2*max.numb.traits.per.SI))){    #columns with significant traits
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

    print(paste("Following parameters being used: common S.I. ",common.SI,"; gene ",gene,"; chr ",chr,"; trait",trait," with peak marker ",QTL,sep = ""))

    #add an object that has all of the names of the NAM founders
    pop.seq <- as.data.frame(as.factor(c("pop01", "pop02", "pop03", "pop04", "pop05", "pop06", "pop07", "pop08", 
                                      "pop09", "pop10", "pop11", "pop12", "pop13", "pop14", "pop15", "pop16", 
                                      "pop18", "pop19", "pop20", "pop21", "pop22", "pop23", "pop24", "pop25", "pop26")))
    founder.names <- as.data.frame(c("B97", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", 
                                 "CML69", "HP301", "IL14H", "KI11", "KI3", "KY21", "M162W", "M37W", "MO18W",                    ### Founder name for Il14H changed to IL14H on  7/13/15
                                 "MS71", "NC350", "NC358", "OH43", "OH7B", "P39", "TX303", "TZI8"))
    NAM.pops <- cbind(pop.seq, founder.names)
    colnames(NAM.pops) <- c("Pop.num", "Pop.Founders")

    if(trait.set == "Carot"){
      effect.estimates.from.entire.trait <- read.table(paste(home.dir,location.of.effect.estimates,"Pop.by.Marker.Effect.Estimates.from.R.",trait,".SI01_2015.txt", sep = ""),head = TRUE)#} 
                                      ####IMPORTANT!!!!!!!!!!!  Requirement for transformed estimates only changed as of June 2015
    }
    if(trait.set == "tocos"){
      effect.estimates.from.entire.trait <- read.table(paste(home.dir,location.of.effect.estimates,"Pop.by.Marker.Effect.Estimates.from.R.",trait,".SI01.txt", sep = ""),head = TRUE)#} 
    }

    markerID <-  substr(effect.estimates.from.entire.trait[,1], 7,10000)
    popID <- substr(effect.estimates.from.entire.trait[,1], 1,5)

    #Obtain transformed effect estimates for this QTL from within master file for this trait
    #transformed.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
                                #as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), ncol(effect.estimates.from.entire.trait)])

    transformed.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
                                as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), 2]))
    colnames(transformed.effect.estimates.from.QTL) <- c("Population","Trans.Est")
    
    lambda = lambda.values[which(lambda.values[,2] == trait),3]
    scalar = if(lambda < 0){-1}else{1}  
    
    print("testing effect signs BEFORE multiplying by -1 if lambda negative.")
    print(transformed.effect.estimates.from.QTL[1:10,2])
    transformed.effect.estimates.from.QTL[,2] = transformed.effect.estimates.from.QTL[,2]*scalar
    print("testing effect signs AFTER multiplying by -1 if lambda negative.")
    print(transformed.effect.estimates.from.QTL[1:10,2])
    transformed.effect.estimates.from.QTL <- merge(transformed.effect.estimates.from.QTL, NAM.pops, by.x = "Population", by.y = "Pop.num")

    #Obtain a vector of all of the GRZM ids of the genes within the interval
    condition.1 <- absolute.final.data.set.FPKM[,3] == paste("chr", chr,sep = "")
    condition.2 <- absolute.final.data.set.FPKM[,4] >= QTL.start.bp
    condition.3 <- absolute.final.data.set.FPKM[,5] <= QTL.stop.bp

    gene.candidate.list.from.SI <- absolute.final.data.set.FPKM[which(condition.1 & condition.2 & condition.3), c(1,2,4,5)]     #used later for ordering effect est data                 

    gene.candidate.list.from.SI <- gene.candidate.list.from.SI[order(gene.candidate.list.from.SI[,3]),]
    list.of.gene.GRZMs.in.support.interval <- as.character(gene.candidate.list.from.SI[,1])                                                      
    list.of.gene.names <- as.character(gene.candidate.list.from.SI[,2])                                                                           


    #########################################################################
    #Run get me my SNPs to obtain the genotypes of these SNPs, read them in below
    if(trait.set == "Carot"){
      imputed.file.location <- paste("QTL_",common.SI,"_imputed_matrix_for_",gene,"\\", sep = "")
      setwd(paste(home.dir,dir.for.5.1.matrices,imputed.file.location,sep = ""))
      file.name <- paste("imputed_RMIP_SNPs_",gene,"_interval_",trait,"_QTL",common.SI,"_chr",chr,".txt",sep = "") 
    }

    if(trait.set == "tocos"){
      imputed.file.location <- paste("QTL_",common.SI,"_imputed_matrix_for_chr",chr,"/", sep = "")
      setwd(paste(home.dir,dir.for.5.1.matrices,imputed.file.location,sep = ""))
      file.name <- paste("imputed_RMIP_SNPs_interval_",trait,"_QTL",common.SI,"_chr",chr,".txt",sep = "") 
    }

    ###CHD note (20150715): File size 311 checks for files that only have header line--will pipe those QTL-trait pairs through only the correlations not involving GWAS SNPs
    if(file.exists(file.name) & file.info(file.name)$size > 311){
      real.file <- read.delim(file.name)
      setwd(paste(home.dir,dir.for.5.1.matrices,imputed.file.location,sep = ""))
      if(trait.set == "Carot"){
        GWAS.SNPs.on.Chr <- GAPIT.Fragment(file.G=paste("imputed_RMIP_SNPs_",gene,"_interval_",trait,"_QTL",common.SI,"_chr",sep = ""),                       
                      file.Ext.G= "txt", file = chr, frag = 1, genoFormat="hapmap", file.fragment=300, Major.allele.zero = TRUE)
      }
      if(trait.set == "tocos"){
        GWAS.SNPs.on.Chr <- GAPIT.Fragment(file.G=paste("imputed_RMIP_SNPs_interval_",trait,"_QTL",common.SI,"_chr",sep = ""),                       
                                        file.Ext.G= "txt", file = chr, frag = 1, genoFormat="hapmap", file.fragment=300, Major.allele.zero = TRUE)
      }
      
      numeric.GWAS.SNPs.data.temp.with.Head <- as.data.frame(GWAS.SNPs.on.Chr$GD)
      for(SNP_pos in unique(numeric.GWAS.SNPs.data.temp.with.Head[,4])){
        row.Ind_this.SNP = which(numeric.GWAS.SNPs.data.temp.with.Head[,4] == SNP_pos)
        NumbOcc = length(row.Ind_this.SNP)
        if(NumbOcc==1){next}
        else if(NumbOcc==2){
          GD.alleles = numeric.GWAS.SNPs.data.temp.with.Head[row.Ind_this.SNP,2]
          print(GD.alleles)
          corresp.GWAS.results = read.table(paste(home.dir,location.of.GWAS.results,trait,"_iterations/Complete_Results_",trait,".txt", sep=""),stringsAsFactors = FALSE,head=TRUE)
          Index.GWAS.hit.allele = which(corresp.GWAS.results[,2]==SNP_pos & as.numeric(corresp.GWAS.results[,4])>4)
          print(Index.GWAS.hit.allele)
          GWAS.hit.allele = corresp.GWAS.results[Index.GWAS.hit.allele,3]
          print(GWAS.hit.allele)
          row.Ind_to.keep = row.Ind_this.SNP[which(GD.alleles == GWAS.hit.allele)]
          if(length(row.Ind_to.keep)==1){
            row.Ind_to.elim = setdiff(row.Ind_this.SNP,row.Ind_to.keep)
            allele.not.hit = numeric.GWAS.SNPs.data.temp.with.Head[row.Ind_to.elim ,2]
            numeric.GWAS.SNPs.data.temp.with.Head = numeric.GWAS.SNPs.data.temp.with.Head[-row.Ind_to.elim,]
            print(paste("Duplicate SNP removed: SNP ",SNP_pos," allele ",allele.not.hit,sep=''))
          }
          else(print(paste("No matching allele found for duplicate SNP at pos ",SNP_pos,sep='')))
        }
        else(print(paste("SNP at pos",SNP_pos,"occurred >2, specifically ",NumbOcc,", times. Please check ",trait,common.SI, sep = '')))
      }
      
      numeric.GWAS.SNPs.data.temp <- t(numeric.GWAS.SNPs.data.temp.with.Head[,-c(1:11)])
      colnames(numeric.GWAS.SNPs.data.temp) <- t(numeric.GWAS.SNPs.data.temp.with.Head[,1])
      final.labels.for.output.table <- paste("Chr_",GWAS.SNPs.on.Chr$GD[,3],"_",GWAS.SNPs.on.Chr$GD[,4], "_bp", sep = "")
      setwd(home.dir)
   
      #Carry out all three correlation functions for triangulation; if all FPKMs are 0 tryCatch loops will ensure correlation is NA (implemented by CHD 20150717)
      final.output.table = correlate.GWAS.SNPs.with.JL.Effect.Estims(numeric.GWAS.SNPs.data.temp = numeric.GWAS.SNPs.data.temp, transformed.effect.estimates.from.QTL = transformed.effect.estimates.from.QTL)
      correlate.FPKMs.with.GWAS.SNPs.and.JL.Effect.Estims(numeric.GWAS.SNPs.data.temp = numeric.GWAS.SNPs.data.temp,transformed.effect.estimates.from.QTL = transformed.effect.estimates.from.QTL,final.output.table=final.output.table)
   
    } else if(file.exists(file.name) & file.info(file.name)$size == 311){
      print(paste("No GWAS SNPs in JL S.I. for common SI", q, ", trait  ",trait,". Thus, only completing 1 leg of tri.: correl(FPKM, JL QTL effect estimate).",sep=''))
      #No GWAS SNPs in JL support interval; carry out only correlation(FPKM, JL QTL effect estimate)
      correlate.FPKMs.with.GWAS.SNPs.and.JL.Effect.Estims(numeric.GWAS.SNPs.data.temp = numeric.GWAS.SNPs.data.temp,transformed.effect.estimates.from.QTL = transformed.effect.estimates.from.QTL,final.output.table=rep(NA,ncol(numeric.GWAS.SNPs.data.temp)+1))
      #Catch added CHD 7/30/2015 to make sure fake imputed file not produced (from previous QTL-trait combination; was observed to happen for gT3, QTL 9)
      tryCatch({rm(paste(output.folder, "Pearson_Between_FPKM_JL_GWAS_",trait, "_",QTL,"_IMPUTED_2015.txt", sep = ""))},
               error = function(e){print("Imputed file already does not exist; good.")})
      } else{
      print(paste("No imputed SNP file exists for trait ",trait," common SI ",common.SI," chr ",chr,". Not even header was printed. Please check!",sep = ""))
    }   

    print(paste("Finished with common SI ",q," trait ",trait,sep=''))
  } # end p trait/ marker combo from list.by.QTL
} #end q QTL intervals
sink()
