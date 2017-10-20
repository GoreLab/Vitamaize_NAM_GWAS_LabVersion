rm(list = ls())
################
#Saved in CBK files as: pleiotropy_network_for_Alex_CBK_markup_dean_mods.r


###Required files:
### (1) Merged tab summary and pleiotropy matrix file generated from script 5b, containing common SI tracker column 
        #IMPORTANT CHECK: trait names in r script vectors must match up with trait names in file column or row headers, else matrices will not be defined

###Output generated:
### (1) Pleiotropy graphs for entire QTL network





#read in pleiotropy matrix
#setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env\\(14)Pleiotropy")
#setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/Pleiotropy/")
setwd("C:/Users/chd45/Documents/Data_From_Mac/Documents/Toco_NAM/")
home.dir = getwd()
pleio.dir <- "/Summary_Tables.Figures/Pleiotropy/Matrix/"
pleio.graph.dir = "/Summary_Tables.Figures/Pleiotropy/Graph/"
trait.set = "tocos" #Options: "carots" or "tocos"
totals.included = TRUE #INPUT REQUIRED*****************
FDR_level = '0.01'
cutoff <- 0.1                           #IMPORTANT!! Input required here
data <- read.table(paste(home.dir,pleio.dir,"Pleiotropic.Output.Matrix.for.",trait.set,".SI01.TAB.SUM.MERGE_with.dT3.redone_sign.corrected.txt",sep =""), head = TRUE,stringsAsFactors=FALSE)

#select subset of columns 
#pleio.matrix <- cbind(data[,1], data[,6], data[,11:19])
#pleio.matrix <- cbind(data[,1], data[,6], data[,21:39]) #for both trait sets
#pleio.matrix <- cbind(data[,1], data[,6], data[,9:18])
pleio.matrix <- cbind(data[,1], data[,6], data[,17:26]) #for size of toco tab summary after sign correction 9/19/2017

#CHD commented out--will use original column names
#colnames(pleio.matrix) <- c("Original_Trait", "Trait_Marker", "ACAR", "PHYF", "BCRY", "LUT", "ZEI", "TOTCAR", "BCAR", "ZEA", "THLYC")

colnames(pleio.matrix)[1:2]=c("Original_Trait","Trait_Marker")



######################################################################
#generate matrix of 1, 0 and -1 specifying if correlation is significant at a=0.01, r=0.495  
#specify the order in which you would like the traits to show up in the network, note that the linear vector will be displayed clockwise                       
#desired.column.order <- c("Trait_Marker", "Original_Trait", "ACAR", "ZEI", "LUT", "PHYF", "THLYC", "ZEA", "BCRY", "BCAR", "TOTCAR")  
#no.totcar.column.order <- c("Trait_Marker", "Original_Trait", "ACAR", "ZEI", "LUT", "PHYF", "THLYC", "ZEA", "BCRY", "BCAR")

#desired.column.order <- c("Original_Trait","Trait_Marker","ACAR_RUV", "ZEI_RUV", "LUT_RUV", "PHYF_RUV", "THLYC_RUV", "ZEA_RUV", "BCRY_RUV", "BCAR_RUV", "TOTCAR_RUV","gT","gT3","aT","aT3","dT","dT3","PC8","totalT","totalT3","totalTocochrs") 
desired.column.order <- c("Original_Trait","Trait_Marker","PC8","dT","gT","aT","totalT","totalTocochrs","totalT3","aT3","gT3","dT3") 
#no.sum.traits.column.order <- c("Original_Trait","Trait_Marker","ACAR_RUV", "ZEI_RUV", "LUT_RUV", "PHYF_RUV", "THLYC_RUV", "ZEA_RUV", "BCRY_RUV", "BCAR_RUV","gT","gT3","aT","aT3","dT","dT3","PC8")
no.sum.traits.column.order <- c("Original_Trait","Trait_Marker","PC8","dT","gT","aT","aT3","gT3","dT3") 

#Reorder the trait names in your data file according to desired.column.order
pleio.matrix <- pleio.matrix[,desired.column.order]    #traits put in desired order for network output, according to Dean
pleio.matrix.no.sum.traits = pleio.matrix[,no.sum.traits.column.order]

#removing the "total carotenoid" trait to make the colnames vector look like no.totcar.column.order in line 13
#pleio.matrix.only.traits <- pleio.matrix[,1:10]
#pleio.matrix.no.sum.traits <- pleio.matrix[,1:15] #doesn't work for new desired column order

#you can use this conditional statement if you want to have flexibility to compare trait vector with and without total traits
if (totals.included == TRUE) {pleio.matrix <- pleio.matrix
  }else{pleio.matrix <- pleio.matrix.no.sum.traits}

#generate placeholder matrix to hold sig pos, sig neg, non sig correlations, and identities
new.pleio.matrix <- matrix(NA, nrow(pleio.matrix), ncol(pleio.matrix))                
colnames(new.pleio.matrix) <- colnames(pleio.matrix)

#columns of interest [,3:11]; rows of interest [1:nrow]
traits <- colnames(pleio.matrix[3:ncol(pleio.matrix)])                                
col.of.interest <- c(3:length(pleio.matrix))                                            

#identify significant correlations
for (i in 1:nrow(pleio.matrix)){
  for (j in col.of.interest){
    if (is.na(pleio.matrix[i,j])) {new.pleio.matrix[i,j] = 0
    }else if(pleio.matrix[i,j]==1){new.pleio.matrix[i,j] = NA
    }else if (pleio.matrix[i,j] > 0) {new.pleio.matrix[i,j] = 1
    }else if (pleio.matrix[i,j] < 0) {new.pleio.matrix[i,j] = -1 
    }else{(print(paste("Pleio matrix [i,j] value is neither 1 (for same trait) nor NA (non-signif correl) nor signif pos or neg, rather is ",pleio.matrix[i,j],". Please check.",sep='')))}
  } # end j columns 
} #end i rows of pleio.matrix
pleio.matrix <- as.data.frame(pleio.matrix)
new.pleio.matrix <- as.data.frame(new.pleio.matrix)

#adds original trait identifiers to pos, neg, non-sig matrix--CHD changed to 1 from 2 because 2 contained peak markers, 1 contained orig traits
rename.pleio.matrix <- cbind(pleio.matrix[,1], new.pleio.matrix [,3:ncol(new.pleio.matrix)])  
  
#######################################################################  
#generate frequency matrices 
  ########################POS Matrix#################################
  #recording total number of QTL for a given trait; recording total number of positive correlations between original and tested trait [i,j]
   sig.pos.freq.matrix <- NULL
   for (i in traits){  
        orig.trait <- rename.pleio.matrix[which((rename.pleio.matrix[,1]) == i),] 
        freq.for.orig.trait <- NULL
       for (j in 2:ncol(orig.trait)){      #where j equals tested trait
       
        #sum QTL in each column (CHD note: works because all cells either NA or 1, or -1 which gets counted below)
        total.QTL <- length(which(is.na(orig.trait[,j])))
        freq.pos <- length(which((orig.trait[,j]) == "1"))
        
        #record QTL frequencies for each tested trait
        record.freq <- total.QTL + freq.pos                  
        freq.for.orig.trait <- c(freq.for.orig.trait, record.freq)
       } #end j loop
   
   #generate matrix for positive correlations    
   sig.pos.freq.matrix <- rbind(sig.pos.freq.matrix, freq.for.orig.trait) 
    } #end i loop
   colnames(sig.pos.freq.matrix) <- traits
   rownames(sig.pos.freq.matrix) <- traits
   
   write.table(sig.pos.freq.matrix,paste(home.dir,pleio.dir,"sig.pos.freq.matrix_sign.corr.txt",sep=''),sep="\t")
  
  ########################NEG Matrix#################################
  #recording total number of QTL for a given trait; recording total number of negative correlations between original and tested trait [i,j]
   sig.neg.freq.matrix <- NULL
   for (i in traits){  
        orig.trait <- rename.pleio.matrix[which((rename.pleio.matrix[,1]) == i),] 
        freq.for.orig.trait <- NULL
       for (j in 2:ncol(orig.trait)){      #where j equals tested trait
       
        #sum QTL in each column
        total.QTL <- length(which(is.na(orig.trait[,j])))
        freq.neg <- length(which((orig.trait[,j]) == "-1"))
        
        #record QTL frequencies for each tested trait
        record.freq <- total.QTL + freq.neg                  
        freq.for.orig.trait <- c(freq.for.orig.trait, record.freq)
       } #end j loop
   
   #generate matrix for negative correlations    
   sig.neg.freq.matrix <- rbind(sig.neg.freq.matrix, freq.for.orig.trait) 
    } #end i loop
   colnames(sig.neg.freq.matrix) <- traits
   rownames(sig.neg.freq.matrix) <- traits
   
   write.table(sig.neg.freq.matrix,paste(home.dir,pleio.dir,"sig.neg.freq.matrix_sign.corr.txt",sep=''),sep="\t")
####################################################################


 #Translation to input for Pat code
 shared_inflo1 <- sig.pos.freq.matrix
 shared_inflo2 <- sig.neg.freq.matrix 


###if fewer than 1% of comparisons show pleiotropy at p<.01, change to 0
###scale by % pleiotropic QTL
###add up sig_pos and sig_neg in both directions
###these 3 matrices are now symmetrical
sig_pos=matrix(NA,length(traits),length(traits))    #matrix of NA, row and column length sized to the total number of traits
sig_neg=matrix(NA,length(traits),length(traits))    #matrix of NA, row and column length sized to the total number of traits
pleiotropy_matrix=matrix(NA,length(traits),length(traits))            

for (i in 1:length(traits)){
  for (j in 1:length(traits)){
    pos=(shared_inflo1[i,j]+shared_inflo1[j,i])/(shared_inflo1[i,i]+shared_inflo1[j,j])           #shared_inflo1 file contains actual number of sig pos corr - name of pos corr matrix
    neg=(shared_inflo2[i,j]+shared_inflo2[j,i])/(shared_inflo1[i,i]+shared_inflo1[j,j])           #shared_inflo2 file contains actual number of sig neg corr - name of neg corr matrix
    sig_pos[i,j]=ifelse((pos+neg)>cutoff,pos,0)     #if sum of pos and neg ratios > cutoff (presumably 0.1), then keep pos corr value, otherwise store 0
    sig_neg[i,j]=ifelse((pos+neg)>cutoff,neg,0)     #if sum of pos and neg ratios > cutoff (presumably 0.1), then keep neg corr value, otherwise store 0
    pleiotropy_matrix[i,j]=sig_pos[i,j]+sig_neg[i,j]   #in pleio matrix cell, store addition of positive and negative ratios
  }}


################################## NOTE: new color matrix generated for each positive, negative and combined network graph   ##########################
###color matrix: col=2 if only positive pleiotropy; col=3 if only negative pleiotropy;
###col=0 if no pleiotropy; col=5 if both positive and negative
#positive network only
col_matrix.pos=matrix(NA,length(traits),length(traits))
for (i in 1:length(traits)){
  for (j in 1:length(traits)){		
    #col_matrix[i,j]=ifelse(sig_pos[i,j]>0&sig_neg[i,j]==0,3,ifelse(sig_pos[i,j]==0&sig_neg[i,j]>0,2,ifelse(sig_pos[i,j]==0&sig_neg[i,j]==0,0,6)))
    col_matrix.pos[i,j]=ifelse(sig_pos[i,j]>0,3,0)
  }}

#negative network only
col_matrix.neg=matrix(NA,length(traits),length(traits))
for (i in 1:length(traits)){
  for (j in 1:length(traits)){		
    #col_matrix[i,j]=ifelse(sig_pos[i,j]>0&sig_neg[i,j]==0,3,ifelse(sig_pos[i,j]==0&sig_neg[i,j]>0,2,ifelse(sig_pos[i,j]==0&sig_neg[i,j]==0,0,6)))
    col_matrix.neg[i,j]=ifelse(sig_neg[i,j]>0,2,0)
  }}
  
#combined network - 3, pos, green; 2, neg, red; 6, combined, magenta; 0, non-sig, white
col_matrix.all=matrix(NA,length(traits),length(traits))
for (i in 1:length(traits)){
  for (j in 1:length(traits)){		
    #col_matrix.all[i,j]=ifelse(sig_pos[i,j]>0&sig_neg[i,j]==0,3,ifelse(sig_pos[i,j]==0&sig_neg[i,j]>0,2,ifelse(sig_pos[i,j]==0&sig_neg[i,j]==0,0,6)))
    #col_matrix.all[i,j]=ifelse(sig_pos[i,j]>0&sig_neg[i,j]==0,4,ifelse(sig_pos[i,j]==0&sig_neg[i,j]>0,2,ifelse(sig_pos[i,j]==0&sig_neg[i,j]==0,0,6))) #changed 3(green) to 5(cyan) or 4 (blue)
    col_matrix.all[i,j]=ifelse(sig_pos[i,j]>0&sig_neg[i,j]==0,"steelblue4",ifelse(sig_pos[i,j]==0&sig_neg[i,j]>0,2,ifelse(sig_pos[i,j]==0&sig_neg[i,j]==0,0,"darkmagenta")))
     }}

  ###################################### NOTE: network graph using line widths based on sig_pos, sig_neg or pleiotropy_matrix matrices are individually generated   #####################
###lty matrix: give dotted line to weaker pleiotropy
#lty_matrix=apply(pleiotropy_matrix,c(1,2),function(x){ifelse(x>0.25,"solid",ifelse(x<0.1&x>0,"dotted",ifelse(x==0,"blank","dashed")))})
library(network)
inflo_network=as.network.matrix(pleiotropy_matrix)
png(paste(home.dir,pleio.graph.dir,trait.set,".genomewide.POS.pleio.cutoff_",cutoff*100,"pct_noSumTraits_tocos_FDR",FDR_level,"_w.dT3.redone_sign.corr.png",sep=""),height=8,width=8,units="in",res=300)
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot.network(inflo_network,label=traits,mode="circle",edge.lwd=pleiotropy_matrix*15,usearrows=F,label.cex=1.7,boxed.labels=F,edge.col=col_matrix,vertex.col=colors,vertex.cex=1.7,vertex.sides=30)
plot.network(inflo_network,label=traits,mode="circle",edge.lwd=sig_pos*25,usearrows=F,label.cex=1.7,label.pos = 0,boxed.labels=F,vertex.col="black",edge.col=col_matrix.pos,vertex.cex=1,vertex.sides=10)

dev.off()

###lty matrix: give dotted line to weaker pleiotropy
#lty_matrix=apply(pleiotropy_matrix,c(1,2),function(x){ifelse(x>0.25,"solid",ifelse(x<0.1&x>0,"dotted",ifelse(x==0,"blank","dashed")))})
#library(network)
inflo_network=as.network.matrix(pleiotropy_matrix)

png(paste(home.dir,pleio.graph.dir,trait.set,".genomewide.NEG.pleio.cutoff_",cutoff*100,"pct_noSumTraits_tocos_FDR",FDR_level,"_w.dT3.redone_sign.corr.png",sep=""),height=8,width=8,units="in",res=300)
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot.network(inflo_network,label=traits,mode="circle",edge.lwd=pleiotropy_matrix*15,usearrows=F,label.cex=1.7,boxed.labels=F,edge.col=col_matrix,vertex.col=colors,vertex.cex=1.7,vertex.sides=30)
plot.network(inflo_network,label=traits,mode="circle",edge.lwd=sig_neg*25,usearrows=F,label.cex=1.7,label.pos = 0,boxed.labels=F,vertex.col="black",edge.col=col_matrix.neg,vertex.cex=1,vertex.sides=10)

dev.off()

###lty matrix: give dotted line to weaker pleiotropy
#lty_matrix=apply(pleiotropy_matrix,c(1,2),function(x){ifelse(x>0.25,"solid",ifelse(x<0.1&x>0,"dotted",ifelse(x==0,"blank","dashed")))})
#library(network)
inflo_network=as.network.matrix(pleiotropy_matrix)
png(paste(home.dir,pleio.graph.dir,trait.set,".genomewide.ALL.pleio.cutoff_",cutoff*100,"pct_noSumTraits_tocos_FDR",FDR_level,"_w.dT3.redone_sign.corr.png",sep=""),height=8,width=8,units="in",res=600)
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot.network(inflo_network,label=traits,mode="circle",edge.lwd=pleiotropy_matrix*15,usearrows=F,label.cex=1.7,boxed.labels=F,edge.col=col_matrix,vertex.col=colors,vertex.cex=1.7,vertex.sides=30)
plot.network(inflo_network,label=traits,displaylabels=F,mode="circle",edge.lwd=pleiotropy_matrix*25,usearrows=F,label.cex=1.7,label.pos = 0,boxed.labels=F,vertex.col="black",edge.col=col_matrix.all,vertex.cex=1,vertex.sides=10)

dev.off()

#colnames(pleiotropy_matrix) = desired.column.order[3:12]
#rownames(pleiotropy_matrix) = desired.column.order[3:12]
#write.table(as.matrix(pleiotropy_matrix),paste(home.dir,pleio.dir,"Pleiotropy.matrix_Final_PercentSharedness_FDR",FDR_level,"_w.dT3.redone.txt",sep=''),sep="\t",row.names = TRUE)
