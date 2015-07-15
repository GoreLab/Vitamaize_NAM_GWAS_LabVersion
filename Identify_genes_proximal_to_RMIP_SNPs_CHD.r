rm(list = ls())

 #############################
 ###NOTE - April 15, 2015
 ###The physical location of a particular SNP in the TOTCAR GWAS exceeds the physical location of chromosome 7 from maizesequence.org, and Robin's data set
   # Chr        bp  allele RMIPx100
#39   7 190031176 G/C_hm1       12
###
### Chromosome terminus from RObin's browser:c(301354135,237068873,232140174,241473504,217872852,169174353,176764762,175793759,156750706,150189435)
### In addition, the terminus of each chromosome in RObin's data set is different from maizesequence.org
###How to reconcile the differences in physical distance?



#Required files:
### (1) Organized FPKM data file
### (2) GWAS results "Complete.results" files

user = 'chd'

#Set the working directory
if(user == 'cbk'){
setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
home.dir <- getwd()
output.dir <- "C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\Candidate gene search\\"
FPKM.file.dir <- "\\(15)Correlated Expression\\"
dir.for.GWAS.results <- "C:\\Users\\ceb19\\Dropbox\\Carotenoid_NAM_Manuscript\\(3) GWAS results\\"
}

if(user == 'chd'){
setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
home.dir <- getwd()
proximal.genes.dir <- "/Expression_Analyses/Tri_Summaries/ProximalGenes/"
FPKM.file.dir <- "/Expression_Analyses/inputsFor5-2/"
GWAS.results.dir <- "/GWAS_Analysis/GWAS_25fam_HMPonly_TASSEL3_alpha01_2015_corr/"
}

absolute.final.data.set.FPKM <- read.table(paste(home.dir,FPKM.file.dir,"FPKM.table.by.gene.annotation.complete.founder.matrix.txt", sep = ""), head = TRUE)
maize.gene.set <- absolute.final.data.set.FPKM[,1:5]

#Complete GWAS Results column in which RMIP data are found
RMIP.column <- 4 

#generate sortable chr column
chr.column <- substr(maize.gene.set[,3],4,6)
maize.gene.set[,3] <- chr.column

#sort maize.gene.set by left bound and chr; CHD changed order of these two lines because otherwise incorr sorting for next.left being sorted
maize.gene.set <- maize.gene.set[order(maize.gene.set[,3]),]
maize.gene.set <- maize.gene.set[order(maize.gene.set[,4]),]

#obtain list of chromosome length in bp, to use as final genomic position for each chromosome in loop
chr.list <- seq(1:10)
#length.list <- c(301433382,237893627,232227970,242029974,217928451,169381756,176810253,175347686,157021084,149627545)    #from maizesequence.org AGPv3, taken 4/15/15
  ## note, the final chromosome position from maizesequence.org is smaller than the physical location of the last gene in some chromosome lists - WHY?????
  ## for the moment, use the right bound of the last gene as the last left.of.next.gene entry
#length.by.chr <- cbind(chr.list, length.list)

#generate intervals from which to locate query SNP, creating column [,6] which will be "left_of_next_gene", starting loop by chromosome
maize.gene.set.with.intervals <- NULL
for (j in chr.list) {
  ##generate data sets for each chromosome (chr.maize.gene.set), separating out maize.gene.set by column [,3]
    chr.maize.gene.set <- maize.gene.set[which(maize.gene.set[,3] == j),]

  ##generate new vector which will indicate the end of interval to be scanned
    left.of.next.gene <- chr.maize.gene.set[-1,4]   #removing first observation from left vector
    left.of.next.gene <- left.of.next.gene - 1      #remove unit of "1" from every observation in new vector
    rt.bound.of.last.gene <- chr.maize.gene.set[nrow(chr.maize.gene.set),5]
    left.of.next.gene <- c(left.of.next.gene, rt.bound.of.last.gene)    #add final observation to vector
    chr.maize.gene.set <- cbind(chr.maize.gene.set, left.of.next.gene)
  
  ##combine chromosome data sets
    maize.gene.set.with.intervals <- rbind(maize.gene.set.with.intervals, chr.maize.gene.set)
} #end j loop of chr in maize.gene.set

###NOTE: maize.gene.set.with.intervals will be smaller than maize.gene.set due to genes on "UNKNOWN" chromosome not being used in interval set



# trait.list <- c("ACAR_RUV", "BCAR_RUV", "BCRY_RUV", "LUT_RUV", "PHYF_RUV", "THLYC_RUV", "TOTCAR_RUV", "ZEA_RUV", "ZEI_RUV")

#trait.list <- c("ACAR_RUV", "BCAR_RUV", "BCRY_RUV", "LUT_RUV", "PHYF_RUV", "THLYC_RUV", "TOTCAR_RUV",  "ZEA_RUV", "ZEI_RUV")

trait.list =  c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs")

#Read in GWAS SNP lists by trait and obtain only those with RMIP > 4
for (i in trait.list){
  setwd(paste(home.dir,GWAS.results.dir,i,"_iterations\\",sep = ""))

    GWAS.results <- read.delim(paste("Complete_Results_",i,".txt", sep=""))

    #obtain only sig RMIP SNPs
    GWAS.results <- GWAS.results[-which(GWAS.results[,RMIP.column] < 5),]
    
    #sort GWAS.results by location [,2] and chr [,1]
    GWAS.results <- GWAS.results[order(GWAS.results[,2]),]
    GWAS.results <- GWAS.results[order(GWAS.results[,1]),]
    
    #generate GWAS SNP query loop to identify hits from maize.gene.set.with.intervals
    all.candidates.for.trait.GWAS.SNPs <- NULL
    for (k in 1:nrow(GWAS.results)){
    
      #specify SNP location of interest
      query.chr <- GWAS.results[k,1]
      query.SNP <- GWAS.results[k,2]
      query.RMIP <- GWAS.results[k,4]
      
      #scan maize.gene.set.with.intervals with chromosome selected and interval bounded by [,4] and [,6]
      gene.hit <- maize.gene.set.with.intervals[which((query.chr == maize.gene.set.with.intervals[,3])&(query.SNP > maize.gene.set.with.intervals[,4])&(query.SNP < maize.gene.set.with.intervals[,5])),]
      
      if(nrow(gene.hit)==0){
        #means that no gene was hit (RMIP SNP was between 2 genes)
        gene.hit = rep(NA,6)
        
        #return gene immediately adjacent to SNP on left
        gene.minus.one <- maize.gene.set.with.intervals[(which((query.chr == maize.gene.set.with.intervals[,3])&(query.SNP > maize.gene.set.with.intervals[,5])&(query.SNP < maize.gene.set.with.intervals[,6]))),]
        
        #return gene immediately adjacent to SNP on right
        gene.plus.one <- maize.gene.set.with.intervals[(which((query.chr == maize.gene.set.with.intervals[,3])&(query.SNP > maize.gene.set.with.intervals[,5])&(query.SNP < maize.gene.set.with.intervals[,6])))+1,]
        
        candidate.class <- c("left.candidate", "exact.hit", "right.candidate")
      }else{
        #assuming we only want to look at "hit" genes, if one was hit (rather than neighbors)
        gene.minus.one <- rep(NA,6)
        gene.plus.one = rep(NA,6)
        #return preceding gene
        gene.minus.one <- maize.gene.set.with.intervals[(which((query.chr == maize.gene.set.with.intervals[,3])&(query.SNP > maize.gene.set.with.intervals[,4])&(query.SNP < maize.gene.set.with.intervals[,5])))-1,]
        gene.minus.one = gene.minus.one[which(!(gene.minus.one[,1] %in% gene.hit[,1])),]
        
        #return following gene
        gene.plus.one <- maize.gene.set.with.intervals[(which((query.chr == maize.gene.set.with.intervals[,3])&(query.SNP > maize.gene.set.with.intervals[,4])&(query.SNP < maize.gene.set.with.intervals[,5])))+1,]
        gene.plus.one = gene.plus.one[which(!(gene.plus.one[,1] %in% gene.hit[,1])),]
        numb.hits = rep("exact.hit",nrow(gene.hit))
        candidate.class <- c("left.candidate", unlist(numb.hits), "right.candidate")
      }
      
      #generate list for a single GWAS SNP
      candidate.list <- rbind(gene.minus.one, gene.hit, gene.plus.one)
      colnames(candidate.list) = c("gene_locus","description","chr","left","right","left.of.next.gene")
        candidate.list <- cbind(candidate.class, candidate.list)
      SNP.info <- c(i, query.chr, query.SNP, query.RMIP)
        #SNP.info.matrix <- as.matrix(rbind(SNP.info, SNP.info, SNP.info))
      numb.rows.this.SNP = length(candidate.class)
      SNP.info.matrix = matrix(rep(SNP.info,numb.rows.this.SNP),ncol = 4,byrow=TRUE)
        colnames(SNP.info.matrix) <- c("trait", "chr", "SNP.position", "RMIP")
        candidate.list.complete.entry <- cbind(SNP.info.matrix, candidate.list)
        
      #generate list for all SNPs in trait list
      all.candidates.for.trait.GWAS.SNPs <- rbind(all.candidates.for.trait.GWAS.SNPs, candidate.list.complete.entry)
      
      } #end k loop of GWAS results
      
      all.candidates.for.trait.GWAS.SNPs = unique(all.candidates.for.trait.GWAS.SNPs)
      write.table(all.candidates.for.trait.GWAS.SNPs, paste(home.dir,proximal.genes.dir,"Genes.proximal.to.RMIP.SNPs.for.",i,".txt", sep=""), sep = "\t", row.names = FALSE)
      
    } #end i trait loop