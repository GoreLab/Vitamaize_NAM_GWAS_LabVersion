rm(list = ls())

 #############################
#Required files:
### (1) Organized FPKM data file
### (2) GWAS results "Complete.results" files

user = 'chd_mac'

if(user == 'chd'){
setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
home.dir <- getwd()
proximal.genes.dir <- "/Expression_Analyses/Tri_Summaries/ProximalGenes/"
FPKM.file.dir <- "/Expression_Analyses/inputsFor5-2/"
GWAS.results.dir <- "/GWAS_Analysis/GWAS_25fam_HMPonly_TASSEL3_alpha01_2015_corr/"
}

if(user == 'chd_mac'){
    setwd("/Users/anybody/Documents/Toco_NAM/")
    home.dir <- getwd()
    proximal.genes.dir <- paste(home.dir,"/Expression/ProximalGenes/",sep='')
    FPKM.file.dir <- paste(home.dir,"/for_FPKMxEffect_comparison/",sep='')
    GWAS.results.dir <- paste(home.dir,"/GWAS/Summaries_by_Trait/",sep='')
}

absolute.final.data.set.FPKM <- read.table(paste(FPKM.file.dir,"FPKM.table.by.gene.ann.complete.matrix.FPKMthreshold.1_filter.by.kernel_across_all.samples.log2trans.txt", sep = ""), head = TRUE)
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
trait.list =  c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs")

#Read in GWAS SNP lists by trait and obtain only those with RMIP > 4
for (i in trait.list){

if(i == "dT3"){
    GWAS.results <- read.table(paste(GWAS.results.dir,"Complete_Results_",i,"_Redone.txt", sep = ""),head = TRUE)#}
}else{
    GWAS.results <- read.table(paste(GWAS.results.dir,"Complete_Results_",i,".txt", sep = ""),head = TRUE)#}
}

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
      write.table(all.candidates.for.trait.GWAS.SNPs, paste(proximal.genes.dir,"Genes.proximal.to.RMIP.SNPs.for.",i,".txt", sep=""), sep = "\t", row.names = FALSE)
      
    } #end i trait loop