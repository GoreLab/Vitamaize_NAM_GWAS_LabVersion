rm(list = ls())

###Required Files:
### (1) Raw P values files from script (5-3)
### (2) FDR adjusted p values by QTL (script 5-3)
### (3) Tabular summary
### 


library(gplots)
library(gridExtra)
trait.set = "tocos" #Options are "Carot" or "tocos"

if(trait.set == "Carot"){ 
#Set the working directory
#setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis\\")
setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\JL_GWAS_FPKM_overlap_analysis_2015_union\\")
home.dir <- getwd()
correlation.output.dir <- "\\Significance Threshold\\FDR_Corrected_Corr_Results\\"
location.of.raw.P.value.results <- correlation.output.dir    #also in specific trait folder
location.of.FDR.P.value.results = correlation.output.dir
proximal.genes.dir <- "C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(15)Correlated Expression\\Candidate gene search\\"
pval.by.QTL.dir = "\\Significance Threshold\\#Pvalue_by_QTL\\"
heatmap.dir <- "\\Significance Threshold\\#Heatmap_by_QTL\\"
tabSummary.path = "C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(31) Tabular Summary Info for 2015 analysis\\LOD scores\\"
}

if(trait.set == "tocos"){
setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
home.dir <- getwd()
correlation.output.dir <- "/Expression_Analyses/GMMS3.1Results_tocos/corrsFrom5.3/"
location.of.raw.P.value.results <- correlation.output.dir    #also in specific trait folder
location.of.FDR.P.value.results = correlation.output.dir
tabSummary.path = "/Tabular_Summaries/"
proximal.genes.dir <- "/Expression_Analyses/Tri_Summaries/ProximalGenes/"
pval.by.QTL.dir = "/Expression_Analyses/GMMS3.1Results_tocos/pval.byQTL/"
heatmap.dir <- "/Expression_Analyses/GMMS3.1Results_tocos/Heatmaps/"
cand.gene.path = "C:/Users/chd45/Documents/Projects/NAM_GP/Inputs/JL/validate_CBK.AEL/"
}

#Read in the tabular summary. Please note that this done only to obtain a list of trait names
#tabular.summary <- read.table(paste(home.dir,"\\(16)Generating Robust Files for Group Review\\Adding family info to JL tabulated results_SI01\\Tab_Sum_of_JL_Carots_Results_for_all traits_SI01_compiled.txt", sep = ""), head = TRUE)
#CHD note 6/22: since only for trait names CHD is using differently formatted file; checked that correct list of traits obtained.
tabular.summary <- read.table(paste(home.dir,tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions_LODscores.txt", sep = ""), head = TRUE)

trait.list <- as.character(unique(tabular.summary[,1]))

#CHD commented out on 6/22--now using trait.list to be consistent with other scripts
#list.of.traits <- c("ACAR_RUV", "BCAR_RUV", "BCRY_RUV", "LUT_RUV", "PHYF_RUV", "THLYC_RUV", "TOTCAR_RUV", "ZEA_RUV", "ZEI_RUV")

cand.gene.matrix = read.table(paste(cand.gene.path,"Tocochromanol_Candidate_Gene_List_GRZMs_R.formatted.txt",sep=""))

top.five.or.something.else <- "5" #Options any integer in quotes

#For loop through the traits
for(trait in trait.list){
  if(trait=='Total_Tocopherols'){trait='totalT'}
  data.proximal.genes = read.table(paste(home.dir,proximal.genes.dir,"Genes.proximal.to.RMIP.SNPs.for.",trait,".txt", sep=""))
    
  setwd(paste(home.dir,location.of.raw.P.value.results, trait,"\\", sep = ""))
  
  #Read in the P-values
  the.raw.P.values <- read.table(paste("Raw.P.values.for.",trait,".txt",sep = ""), head = TRUE)
  the.FDR.P.values =  read.table(paste("FDR.Adjusted.P.values.for.",trait,"_by_Q.txt",sep = ""), head = TRUE)
  
  #Read in the correlations
  correlations = read.table(paste("Correlations.for.",trait,".txt",sep=""),head=TRUE)
  correlations = as.matrix(correlations)

  #create subdirectories  #note CBK moved this out of QTL loop below on 7/5/15 to allow directory creation to occur once
  trait.subdirectory <- paste(trait, "/", sep = "")
  dir.create(paste(home.dir, pval.by.QTL.dir, trait.subdirectory, sep = ""))
  pval.trait.location <-(paste(home.dir, pval.by.QTL.dir, trait.subdirectory, sep = ""))
  
  
  trait.subdirectory <- paste("P-value_Heatmaps_for_", trait, "\\", sep = "")
  dir.create(paste(home.dir, heatmap.dir, trait.subdirectory, sep = ""))
  heatmap.trait.location <-(paste(home.dir, heatmap.dir, trait.subdirectory, sep = ""))
  
  #setwd(paste(home.dir, heatmap.dir, sep = ""))         #commented out by CBK on 7/5/15, setting up loop in wrong directory
  #For loop through each QTL
  for(QTL in as.character(unique(the.FDR.P.values[,ncol(the.FDR.P.values)]))){
    the.raw.P.values.for.this.QTL <- the.raw.P.values[which(the.raw.P.values[,ncol(the.raw.P.values)] == QTL),]
    the.FDR.P.values.for.this.QTL <- the.FDR.P.values[which(the.FDR.P.values[,ncol(the.FDR.P.values)] == QTL),]
    
    common.SI <- tabular.summary[which((as.character(tabular.summary[,1]) == trait) & (as.character(tabular.summary[,6]) == QTL)),12]   #added by CBK on 7/5/15 to facilitate data assembly 

    #write.table(the.FDR.P.values.for.this.QTL,paste("FDR.P.values.for.",trait,"_QTL.with.peak.at.",QTL,".txt",sep=''),sep='\t',row.names=FALSE)   #commented out by CBK on 7/5/15
    write.table(the.FDR.P.values.for.this.QTL,paste(pval.trait.location, "FDR.P.values.for.",trait,".QTL",common.SI,".QTL.",QTL,".txt",sep=''),sep='\t',row.names=FALSE)
        
    master.top.P.value.list <- NULL

    #For loop through each time point
    for(time.point in 2:7){
      #Sort from smallest to largest P-value
      the.FDR.P.values.for.this.QTL <- the.FDR.P.values.for.this.QTL[order(the.FDR.P.values.for.this.QTL[,time.point]),]
      
      the.top.P.values.for.this.QTL <- the.FDR.P.values.for.this.QTL[1:as.numeric(top.five.or.something.else),]
      
      #generate merged file #added by CBK on 7/5/15
      the.top.P.values.for.this.QTL.mod <- cbind(the.top.P.values.for.this.QTL, rep(colnames(the.top.P.values.for.this.QTL)[time.point], nrow(the.top.P.values.for.this.QTL)))
      master.top.P.value.list <- rbind(master.top.P.value.list, the.top.P.values.for.this.QTL.mod)
      

      #Write a text file that will give information about the
      #write.table(the.top.P.values.for.this.QTL, paste("FPKM.matrix.for.top.",top.five.or.something.else,".Genes.Associated.with.",trait,".",QTL,".on.", the.raw.P.values.for.this.QTL[1,(ncol(the.raw.P.values.for.this.QTL)-3)],".",
      #                               substr(colnames(the.top.P.values.for.this.QTL)[time.point], start = 2, stop = 1000), ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
      write.table(the.top.P.values.for.this.QTL, paste(pval.trait.location, "FPKM.matrix.for.top.",top.five.or.something.else,".Genes.Assoc.with.",trait,".QTL",common.SI,".QTL.",QTL, ".",colnames(the.top.P.values.for.this.QTL)[time.point],".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
      
      #Start a pdf device
     #pdf(paste("Top.",top.five.or.something.else,".Genes.Associated.with.",trait,".",QTL,".on.", the.raw.P.values.for.this.QTL[1,(ncol(the.raw.P.values.for.this.QTL)-3)],".pdf" , sep = ""), width = 14)      #commented out by CBK on 7/5/2015 to add in common SI info
      pdf(paste(heatmap.trait.location, "Top.",top.five.or.something.else,".Genes.for.",trait,".QTL",common.SI,".QTL.",QTL,".",colnames(the.top.P.values.for.this.QTL)[time.point],".pdf" , sep = ""), width = 14)       #added by CBK on 7/5/2015 to add in common SI info

      #Make a heatmap of the "top 5" genes with the strongest association
      data.for.heatmap <- as.matrix(the.top.P.values.for.this.QTL[,2:7])
      data.for.heatmap[which(is.na(data.for.heatmap))] = 1
      #Take the -log10 P-values so that the stronger correlations are darker in the heatmap
      data.for.heatmap <- -log10(data.for.heatmap)
      rownames(data.for.heatmap) <- the.top.P.values.for.this.QTL[,1]
      
      abs.max <- max(abs(data.for.heatmap))
      increments <- (2*abs.max)/20
      breaks = seq(-abs.max, abs.max, by = increments)  
     if(increments == 0){
       print(paste("By increments criterion, no non-zero log10 pvalues for this timepoint ",time.point,trait,common.SI,QTL,sep=''))
       next
     }

      heatmap.2(data.for.heatmap, dendrogram = "none", Rowv = FALSE, Colv = "Rowv", trace = "none", 
                key = TRUE, col = bluered, ylab = NULL, breaks = breaks,  
                margins = c(5,10), cexRow = 0.7, cexCol = 0.7, main = paste("Top ", top.five.or.something.else, " Correlated Genes at ", 
                                                                            substr(colnames(the.top.P.values.for.this.QTL)[time.point], start = 2, stop = 1000), "(-log10 P-value)" , sep = ""))      
       dev.off()   
    }#End for loop through each time point

    correlations.for.this.QTL = correlations[which(correlations[,ncol(correlations)] == QTL),]
    top.correls.this.QTL = as.data.frame(matrix(NA,nrow=nrow(master.top.P.value.list),ncol=ncol(correlations.for.this.QTL)))
    top.raw.P.this.QTL = as.data.frame(matrix(NA,nrow=nrow(master.top.P.value.list),ncol=ncol(the.raw.P.values.for.this.QTL)))
    for(row in 1:nrow(top.correls.this.QTL)){
      top.gene = master.top.P.value.list[row,1]
      if(is.na(top.gene)){
        top.correls.this.QTL[row,] = rep(NA,ncol(top.correls.this.QTL))
        top.raw.P.this.QTL[row,] = rep(NA,ncol(top.raw.P.this.QTL))
      }else{    
      for(col in 1:ncol(top.correls.this.QTL)){
        top.correls.this.QTL[row,col] = correlations.for.this.QTL[which(correlations.for.this.QTL[,1] == top.gene),col]
      }
      for(col2 in 1:ncol(top.raw.P.this.QTL)){
        top.raw.P.this.QTL[row,col2] = the.raw.P.values.for.this.QTL[which(the.raw.P.values.for.this.QTL[,1] == top.gene),col2]
      }#end for col loop 
    }#end else
    }#end for row loop
    colnames(top.correls.this.QTL)=colnames(correlations.for.this.QTL)
    colnames(top.raw.P.this.QTL)=colnames(the.raw.P.values.for.this.QTL)
    
    #print out merged data table of all ranked timepoints
     colnames(master.top.P.value.list)[15] <- "timepoint.used.for.ranking"
    #CHD changed row names to false on 7/7, otherwise col names shifted by one
     write.table(master.top.P.value.list, paste(pval.trait.location, "FPKM.matrix.for.top.",top.five.or.something.else,".Genes.Assoc.with.",trait,".QTL",common.SI,".QTL.",QTL, ".all.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
      
     #Start a pdf device for the merged data file ### CBK added on 7/5/15
     #pdf(paste("Top.",top.five.or.something.else,".Genes.Associated.with.",trait,".",QTL,".on.", the.raw.P.values.for.this.QTL[1,(ncol(the.raw.P.values.for.this.QTL)-3)],".pdf" , sep = ""), width = 14)      #commented out by CBK on 7/5/2015 to add in common SI info
      pdf(paste(heatmap.trait.location, "Top.",top.five.or.something.else,".Genes.for.",trait,".QTL",common.SI,".QTL",QTL,"all.pdf" , sep = ""), width = 14)       #added by CBK on 7/5/2015 to add in common SI info

      #Added by CHD 7/7: print a tabular summary of genes appearing in which time points, and whether RMIP SNP is within or adjacent to gene
    allTimes.Tab.summ = matrix(NA, nrow = length(unique(master.top.P.value.list[,1])), ncol=26)
    colnames(allTimes.Tab.summ) = c("GRZM_ID","annot","RGv2_start","RGv2_stop","cand?","hit.by.RMIP.SNP","adjac.to.RMIP.SNP","kbp.from.adjac.SNP","r.12_DAP","Fp.12_DAP","Rp.12_DAP","r.16_DAP","Fp.16_DAP","Rp.16_DAP","r.20_DAP","Fp.20_DAP","Rp.20_DAP","r.24_DAP","Fp.24_DAP","Rp.24_DAP","r.30_DAP","Fp.30_DAP","Rp.30_DAP","r.36_DAP","Fp.36_DAP","Rp.36_DAP")
    for(i in (1:nrow(allTimes.Tab.summ))){
      thisGene = as.character(unique(master.top.P.value.list[,1])[i])
      if(is.na(thisGene)){next}
      #insert GRZM ID
      allTimes.Tab.summ[i,"GRZM_ID"]= thisGene
      #insert functional annotation
      allTimes.Tab.summ[i,"annot"]= as.character(unique(master.top.P.value.list[which(master.top.P.value.list[,1]==thisGene),10]))
      #insert gene start pos
      allTimes.Tab.summ[i,"RGv2_start"]= unique(master.top.P.value.list[which(master.top.P.value.list[,1]==thisGene),12])
      #insert gene stop pos
      allTimes.Tab.summ[i,"RGv2_stop"]= unique(master.top.P.value.list[which(master.top.P.value.list[,1]==thisGene),13])

      #indicate candidate gene or not
      if(thisGene %in% cand.gene.matrix[,2]){
        allTimes.Tab.summ[i,"cand?"] = "Yes"
      }else{allTimes.Tab.summ[i,"cand?"] = "."}

      FDR.p_timepoints.thisGene = master.top.P.value.list[which(master.top.P.value.list[,1]==thisGene),]
      corr_timepoints.thisGene = top.correls.this.QTL[which(top.correls.this.QTL[,1]==thisGene),]
      raw.P.thisGene = top.raw.P.this.QTL[which(top.raw.P.this.QTL[,1]==thisGene),]
      for(j in (1:nrow(FDR.p_timepoints.thisGene))){
        FDR.p_thisTimepoint = FDR.p_timepoints.thisGene[j,]
        corr_thisTimepoint = corr_timepoints.thisGene[j,]
        raw.P_thisTimepoint = raw.P.thisGene[j,]
        timeRanked = as.character(FDR.p_thisTimepoint[,15])
        pval.at.timeRanked = FDR.p_thisTimepoint[,timeRanked]
        corr.at.timeRanked = corr_thisTimepoint[,timeRanked]
        raw.P.at.timeRanked = raw.P_thisTimepoint[,timeRanked]
        timeRanked = substr(timeRanked,start=2,stop=10)
        geneIndex = which(allTimes.Tab.summ[,1] == thisGene,)
        Fp.timeRanked.colname = paste("Fp.",timeRanked,sep='')
        r.timeRanked.colname = paste("r.",timeRanked,sep='')
        Rp.timeRanked.colname = paste("Rp.",timeRanked,sep='')
        allTimes.Tab.summ[geneIndex,Fp.timeRanked.colname]= round(as.numeric(pval.at.timeRanked), digits=3)
        allTimes.Tab.summ[geneIndex,r.timeRanked.colname] = round(as.numeric(corr.at.timeRanked),digits=3)
        allTimes.Tab.summ[geneIndex,Rp.timeRanked.colname] = round(as.numeric(raw.P.at.timeRanked),digits=3)
      }
      
      if(thisGene %in% data.proximal.genes[,6]){
        hits.thisGene = data.proximal.genes[which(data.proximal.genes[,6]==thisGene),]
        for(k in (1:nrow(hits.thisGene))){
          thisHit = hits.thisGene[k,]
          if(thisHit[5]=="exact.hit"){allTimes.Tab.summ[geneIndex,"hit.by.RMIP.SNP"]=as.character(thisHit[,3])}
          if(thisHit[5]=="left.candidate"){
            allTimes.Tab.summ[geneIndex,"adjac.to.RMIP.SNP"]=as.character(thisHit[,3])
            allTimes.Tab.summ[geneIndex,"kbp.from.adjac.SNP"]=round(as.numeric((as.numeric(thisHit[10])-as.numeric(thisHit[3]))/1000),digits=3)           
          }
          if(thisHit[5]=="right.candidate"){
            allTimes.Tab.summ[geneIndex,"adjac.to.RMIP.SNP"]=as.character(thisHit[,3])
            allTimes.Tab.summ[geneIndex,"kbp.from.adjac.SNP"]=round(as.numeric((as.numeric(thisHit[9])-as.numeric(thisHit[3]))/1000),digits=3)
        }
      }
    }
    }
    the.FDR.P.values.for.this.QTL[order(the.FDR.P.values.for.this.QTL[,time.point]),]
    allTimes.Tab.summ.ordered = allTimes.Tab.summ[order(allTimes.Tab.summ[,"RGv2_start"]),]
    allTimes.Tab.summ.ordered[is.na(allTimes.Tab.summ.ordered)] = "."
    grid.table(allTimes.Tab.summ.ordered, gpar.coretext = gpar(fontsize=4), gpar.coltext = gpar(fontsize=3), padding.h=unit(1, "mm"), padding.v=unit(1, "mm"), show.rownames = FALSE,
               h.even.alpha = 1, h.odd.alpha = 1, v.even.alpha = 0, v.odd.alpha = 0.4, gpar.corefill = gpar(fill = 'blue',col = 'black'))
    
    
    #Make a heatmap of the "top 5" genes with the strongest association
      data.for.heatmap <- as.matrix(master.top.P.value.list[,2:7])
      data.for.heatmap[which(is.na(data.for.heatmap))] = 1
      #Take the -log10 P-values so that the stronger correlations are darker in the heatmap
      data.for.heatmap <- -log10(data.for.heatmap)
      rownames(data.for.heatmap) <- master.top.P.value.list[,1]
      
      abs.max <- max(abs(data.for.heatmap))
      increments <- (2*abs.max)/20
    if(increments == 0){
      print(paste("By increments criterion, no timepoints with non-zero log10 pvalues for",trait,common.SI,QTL,sep=''))
      next
    }
      breaks = seq(-abs.max, abs.max, by = increments)  
    
    if(length(unique(c(data.for.heatmap)))==1){
      print(paste("By all values are same criterion, no timepoints with non-zero log10 pvalues for",trait,common.SI,QTL,sep=''))
      next
    }
    
    rowIdx_to_remove = NULL
    for(row in (1:nrow(data.for.heatmap))){
      if(length(unique(c(data.for.heatmap[row,])))==1 && unique(c(data.for.heatmap[row,])==0)){
        rowIdx_to_remove = c(rowIdx_to_remove,row)
      }
    }
    if(!(length(rowIdx_to_remove)==0)){
    data.for.heatmap = data.for.heatmap[-rowIdx_to_remove,]
    }
    
    if((nrow(data.for.heatmap)==0 || is.null(nrow(data.for.heatmap)))){
      print(paste("By nrow & unique criteria after removing all-0 rows, no timepoints with non-zero log10 pvalues for",trait,common.SI,QTL,sep=''))
      next
    }else{
    #CHD changed main to all time points rather than substr(colnames(the.top.P.values.for.this.QTL)[time.point], start = 2, stop = 1000) (was using 36_DAP, not approp for merged map)
      heatmap.2(data.for.heatmap, dendrogram = "none", Rowv = FALSE, Colv = "Rowv", trace = "none", 
                key = TRUE, col = bluered, ylab = NULL, breaks = breaks,  
                margins = c(5,10), cexRow = 0.7, cexCol = 0.7, main = paste("Top ", top.five.or.something.else, " Correlated Genes at all time points, -log10 FDR(by Q) P-value" , sep = ""))      
    }  
      

    dev.off()
    #Turn off the pdf device
  


  }#End for loop through each QTL



} #End for(trait in trait.list)