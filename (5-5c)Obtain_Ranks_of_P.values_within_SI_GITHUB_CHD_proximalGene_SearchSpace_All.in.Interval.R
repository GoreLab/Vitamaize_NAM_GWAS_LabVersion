###Note to Cathy: dir.for.compil.of.tri.summ now needs to be specified rather than pval.by.QTL.dir
#CHD used dir.for.5.2.ordered (created in 5-2) to output p-values to same directory as imputed and ordered files needed for compilation of tri summaries
#can still use separate dir if desired

rm(list = ls())


GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure.mod <-
  function(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH"){
    #Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
    #Output: PWIP, number.of.significant.SNPs
    #Authors: Alex Lipka and Zhiwu Zhang 
    # Last update: May 5, 2011 
    ##############################################################################################
    #Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest
    if(is.null(PWI)){
      PWIP=NULL
      number.of.significant.SNPs = 0
    }
    if(!is.null(PWI)){  
      if(length(which(is.na(PWI)))>0){
        print(str(PWI))
        print(length(which(is.na(PWI))))
        PWI[is.na(PWI)] = 1
        #PWI[which(is.na(PWI))] = 1   #This line was in Alex's code but does not work for this data set (perhaps due to type)
        print(PWI)
      }
      for(i in 1:ncol(PWI)){
        #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
        #Right column: Simes corrected p-value
        res <- mt.rawp2adjp(as.vector(PWI[,i]), FDR.Procedure)
        
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


#########################################################################################################################################################################################################
###Required Files:
### (1) Raw P values files from script (5-3)

library(gplots)
library(multtest)

user = "cbk"              #Options are "cbk" or "chd" (maintains file path differences to keep script easily adaptable)
pvals.to.output = "search.range" #Options are 
                          #"1gene": directly hit gene + 1 on either side, or if SNP falls between 2 genes only those 2
                          #"search.range": certain physical distance in each direction from RMIP SNP, even if outside of individual or common S.I.
                          #"all.in.interval": all RMIP SNPs falling within individual trait-JL QTL SI within a common S.I. 
search.range <- 100000    #only is used if pvals.to.output == "search.range"

if(user == 'cbk'){
  setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
  home.dir <- getwd()
  correlation.output.dir <- "\\(15)Correlated Expression\\Significance_Threshold_union\\FDR_Corrected_Corr_Results\\"
  location.of.raw.P.value.results <- correlation.output.dir    #also in specific trait folder
  location.of.FDR.P.value.results = correlation.output.dir
  proximal.genes.dir <- "\\(15)Correlated Expression\\Candidate gene search\\"
  dir.for.compil.of.tri.summ = "\\(15)Correlated Expression\\Results_from_R_FPKM1_logtyes\\"
  pval.by.QTL.dir = "\\(15)Correlated Expression\\Significance_Threshold_union\\Pvalue_by_QTL\\"
  trait.set = "Carot" #Options are "Carot" or "tocos"
  tabSummary.path = "\\(31) Tabular Summary Info for 2015 analysis\\"
  location.of.GWAS.results <- "\\(10)GWAS Analysis\\RUV GWAS 25fam_alldata_alpha01_2015_FINAL\\"
  column.of.cSI.in.tab.sum <- 14
  trait.col.common.key <- 1
  chr.col.common.key <- 2
  pos.col.common.key <- 3
  
  #Read in the tabular summary; used below to obtain a list of trait names
  #tabular.summary <- read.table(paste(home.dir,"\\(16)Generating Robust Files for Group Review\\Adding family info to JL tabulated results_SI01\\Tab_Sum_of_JL_Carots_Results_for_all traits_SI01_compiled.txt", sep = ""), head = TRUE)
  tabular.summary <- read.table(paste(home.dir,tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_GWAS_FamPVE_common_SI_recsuppregions_LODscores_20150612.txt", sep = ""), head = TRUE)

}

if(user == 'chd'){
  setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
  home.dir <- getwd()
  correlation.output.dir <- "/Expression_Analyses/GMMS3.1Results_tocos/corrsFrom5.3/"
  location.of.raw.P.value.results <- correlation.output.dir    #also in specific trait folder
  location.of.FDR.P.value.results = correlation.output.dir
  proximal.genes.dir <- "/Expression_Analyses/Tri_Summaries/ProximalGenes/"
  dir.for.compil.of.tri.summ = "/Expression_Analyses/GMMS3.1Results_tocos/orderedFrom5.2/" 
  trait.set = "tocos" #Options are "Carot" or "tocos"
  tabSummary.path = "/Tabular_Summaries/"
  location.of.GWAS.results = "/GWAS_Analysis/GWAS_25fam_HMPonly_TASSEL3_alpha01_2015_corr/"
  fdr.adjusted.for.100kb.dir = correlation.output.dir
  
  #Read in the tabular summary; used below to obtain a list of trait names
  tabular.summary <- read.table(paste(home.dir,tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_SI_with_GWAS_SNPs_common_SI_20150511_recsuppregions_LODscores.txt", sep = ""), head = TRUE)
}

trait.list <- as.character(unique(tabular.summary[,1]))

#CHD added 7/5: Key file so that output files can be labeled by common S.I. number.
#Note CHD changed to use file with suffix dupSNPremoved in which duplicate SNPs across HapMap v.1 and 2 are consolidated (only allele column was different); otherwise proximal genes were being printed in duplicate
<<<<<<< HEAD
RMIP.commonSI.key = read.table(paste(home.dir,tabSummary.path,"Complete_Results_allTraits_RMIPgt4_with_QTLnumber.txt", sep = ""), head = TRUE)
=======
RMIP.commonSI.key = read.table(paste(home.dir,tabSummary.path,"Complete_Results_allTraits_RMIP_test_with_QTLnumber.txt", sep = ""), head = TRUE)

if(user == 'cbk'){
Common.SI.array = read.table(paste(home.dir,dir.for.compil.of.tri.summ,"Common_SI_array_for_tri_auto_June_2015.txt", sep = ""), head = TRUE)       #needed only by cbk
}
>>>>>>> ebcdcbfb5d502380201104b1474ad81ab4064123

#For loop through each common S.I.(39 for carot, 48 for toco)
for(cSI in unique(tabular.summary[,column.of.cSI.in.tab.sum])){
  print(paste("For common support interval number ",cSI,":",sep=''))
  
  #Initialize master p-value objects that will hold results compiled within a common S.I. (compiled across trait/JL marker pairs for all.in.interval, across RMIP SNPs for 1gene or search.space)
  FDR.pval.this.cSI_allTraits = NULL
  fits.criteria.raw.master = NULL
  fits.criteria.FDR.master = NULL
  proximal.summary.master = NULL
  
  #For loop through the traits
  for(trait in trait.list){
  
    if(user == 'chd'){
    setwd(paste(home.dir,location.of.raw.P.value.results, trait,"\\", sep = ""))
    }
    
    if(user == 'cbk'){
     setwd(paste(home.dir,location.of.raw.P.value.results, sep = ""))
     }
     
    #Read in the raw and FDR P-values
    the.raw.P.values <- read.table(paste("Raw.P.values.for.",trait,".txt",sep = ""), head = TRUE, stringsAsFactors = FALSE)
    the.FDR.P.values =  read.table(paste("FDR.Adjusted.P.values.for.",trait,"_by_Q.txt",sep = ""), head = TRUE, stringsAsFactors = FALSE)
    
    #Change dir to master output dir for this common S.I.; final files from this script will be placed here.
<<<<<<< HEAD
    setwd(paste(home.dir,dir.for.compil.of.tri.summ,"QTL_",cSI,"_tri.files_FPKM1_logtyes/", sep=''))
=======
>>>>>>> ebcdcbfb5d502380201104b1474ad81ab4064123
    
    if(user == 'cbk'){
      gene.abbrev <- Common.SI.array[which(Common.SI.array[,1] == cSI),5]                                                    #cbk modified
      setwd(paste(home.dir,dir.for.compil.of.tri.summ,"QTL_",cSI,"_imputed.ordered.tri.files.for.",gene.abbrev, sep=''))     #cbk modified
    }  
      
    if(pvals.to.output == "all.in.interval"){
      
      #Create list of peak markers for this trait
      all_JL.markers.this.trait = as.character(unique(the.FDR.P.values[,ncol(the.FDR.P.values)]))
      
      #Extract tabular summary rows corresponding to this trait-common S.I. combination
      tabular.summary.this.cSI = tabular.summary[which(
          tabular.summary[,1] == trait &
          tabular.summary[,6] %in% all_JL.markers.this.trait & 
          tabular.summary[,column.of.cSI.in.tab.sum] == cSI),]
      
      #If no hit for this trait within this common S.I., go to next trait.
      if(length(unique(tabular.summary.this.cSI[,6]))==0){
        next
      }else { 
        print(paste("Now processing, for all.in.interval, trait ",trait,sep=''))
        
        #Iterate through peak markers for this trait-SI pair (determine marker IDs ("S_{pos}") from tabular summary subset, 'Peak_Marker' column)
        for(JL.marker in unique(tabular.summary.this.cSI[,6])){
          
          #Extract FDR p-values for this marker from master p-value file compiled per trait across common S.I.s
          FDR.pval.this.cSI_1trait = the.FDR.P.values[which(the.FDR.P.values[,ncol(the.FDR.P.values)] == JL.marker),]
          
          #Check common SI number to confirm still in correct interval
          cSI.check = unique(tabular.summary[which(tabular.summary[,6]==JL.marker),column.of.cSI.in.tab.sum])
          if(!(cSI == cSI.check)){
            print(paste("Something went wrong for peak JL marker",JL.marker," for trait", trait,". Results being calculated outside of correct common support interval.",sep=''))
            }
        
        #Add trait column identifier
        FDR.pval.this.cSI_1trait = cbind(FDR.pval.this.cSI_1trait,rep(trait,nrow(FDR.pval.this.cSI_1trait)))
        colnames(FDR.pval.this.cSI_1trait)[ncol(FDR.pval.this.cSI_1trait)] = "trait"
        
        #Append this trait to master file for this common S.I.
        FDR.pval.this.cSI_allTraits = rbind(FDR.pval.this.cSI_allTraits,FDR.pval.this.cSI_1trait)
        
        #####If want to output by single trait:
        #write.table(the.FDR.P.values.for.this.QTL,paste("FDR.P.values_",trait,"_all.in.interval_QTL",cSI,".txt",sep=''),sep='\t',row.names=FALSE)
        
        }#end loop through peak JL markers
        }#end loop indicating non-zero number of peak JL markers to process
    } # end if("all.in.interval")
    
    #Extract RMIP common S.I. key rows corresponding to this trait-common S.I. combination
    RMIP.SNPs.this.trait_cSI = RMIP.commonSI.key[which(RMIP.commonSI.key[,trait.col.common.key]==trait & RMIP.commonSI.key[,6]==cSI),]
    
    #Remove duplicate chr,pos pairs but also combine alleles so keep record of duplicate markers between HapMap v1 and v2
    #aggregate(allele ~ Chr+bp, data = RMIP.SNPs.this.trait_cSI, FUN = cat)
    
    #If no RMIP SNPs for this trait within this common S.I., go to next trait.
    if(nrow(RMIP.SNPs.this.trait_cSI) == 0){
      next
    }else{
      if(pvals.to.output == "search.range"){ #end catch for no JL/GWAS overlap, start 100kb (or other specified physical distance) search
        print(paste("Now processing, for search.range, trait ",trait,sep=''))
        
        count <- 0
        #Iterate through RMIP SNPs for this trait-SI pair
        for(RMIP.SNP in (1:nrow(RMIP.SNPs.this.trait_cSI))){
          #Determine RMIP SNP chr and pos from RMIP common S.I. key subset, columns 1 and 2
          chr = RMIP.SNPs.this.trait_cSI[RMIP.SNP, chr.col.common.key]
          pos = RMIP.SNPs.this.trait_cSI[RMIP.SNP, pos.col.common.key]
          
          #Identify search space for +/-{specified physical distance} region
          lower.search.bound <- as.numeric(pos - search.range)
          upper.search.bound <- as.numeric(pos + search.range)                                     
          
          #Identify genes within the search space - raw
          the.raw.P.values.subset <- the.raw.P.values[which(substr(the.raw.P.values[,11], 4, 5) == chr),]
          condition.1 <- the.raw.P.values.subset[,12] > lower.search.bound 
          condition.2 <- the.raw.P.values.subset[,13] < upper.search.bound
          
          #Extract raw p-values for these genes
          fits.criteria.raw.this.SNP <- the.raw.P.values.subset[which(condition.1 & condition.2),] 
          
          #Add trait column identifier
          fits.criteria.raw.this.SNP = cbind(fits.criteria.raw.this.SNP,rep(trait,nrow(fits.criteria.raw.this.SNP)))
          colnames(fits.criteria.raw.this.SNP)[ncol(fits.criteria.raw.this.SNP)] = "trait"
        
          #Append these values to master file for this common S.I.
          fits.criteria.raw.master = rbind(fits.criteria.raw.master,fits.criteria.raw.this.SNP)
          
          print(paste("----------------And now, we are calculating the FDR-adjusted P-values BY 100kb for ", trait, "------------------------", sep = ""))
          #Run the FDR correction
          #if no genes within 100kb skip this SNP so that GAPIT.BH will not hit error
          if(nrow(fits.criteria.raw.this.SNP)==0){
            print(paste("No genes within 100kb of SNP ",RMIP.SNP, "thus no FDR p-values calculated.",sep=''))
            next
          #if only 1 row in the p-val matrix (due to only 1 gene in search space), capture values 'as is' (i.e. raw) so that GAPIT.BH will not hit error
          }else if(nrow(fits.criteria.raw.this.SNP)==1){
              FDR.by100kb.Adjusted.P.values = fits.criteria.raw.this.SNP[,2:9]
              print(paste("No adjustment made for p-values based on 100kb FDR because only 1 gene within 100kb, thus no multiple testing."))
          }else{
            FDR.by100kb.Adjusted.P.values <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure.mod(PWI = fits.criteria.raw.this.SNP[,2:9])
          }
          
            FDR.by100kb.P.value.results <- cbind(fits.criteria.raw.this.SNP[,1], FDR.by100kb.Adjusted.P.values, 
                                               fits.criteria.raw.this.SNP[,10:14])
            colnames(FDR.by100kb.P.value.results)[1:14] <- colnames(fits.criteria.raw.this.SNP)[1:14]

            #Add trait column identifier
            FDR.by100kb.P.value.results = cbind(FDR.by100kb.P.value.results,rep(trait,nrow(FDR.by100kb.P.value.results)))
            colnames(FDR.by100kb.P.value.results)[ncol(FDR.by100kb.P.value.results)] = "trait"

          #Append the support interval-wise FDR corrected P-values
          if(count == 0){
            FDR.adjusted.P.values.all.SNPs.this.QTL <- FDR.by100kb.P.value.results
          }else{
            FDR.adjusted.P.values.all.SNPs.this.QTL <- rbind(FDR.adjusted.P.values.all.SNPs.this.QTL,FDR.by100kb.P.value.results) 
          }
          count <- count + 1
          
          #####If want to use FDR p-values rather than raw:
          #Identify genes within the search space - FDR
          #the.FDR.P.values.subset <- the.FDR.P.values[which(substr(the.FDR.P.values[,11], 4, 5) == chr),]
          #condition.1 <- the.FDR.P.values.subset[,12] > lower.search.bound 
          #condition.2 <- the.FDR.P.values.subset[,13] < upper.search.bound
          
          #Extract raw p-values for these genes
          #fits.criteria.FDR.this.SNP <- the.FDR.P.values.subset[which(condition.1 & condition.2),] 
          
          #Add trait column identifier
          #fits.criteria.FDR.this.SNP = cbind(fits.criteria.FDR.this.SNP,rep(trait,nrow(fits.criteria.FDR.this.SNP)))
          #colnames(fits.criteria.FDR.this.SNP)[ncol(fits.criteria.FDR.this.SNP)] = "trait"
          
          #Append these values to master file for this common S.I.
          #fits.criteria.FDR.master = rbind(fits.criteria.FDR.master,fits.criteria.FDR.this.SNP)
          
          } #end loop through RMIP SNPs
        
        #####If want to output by single trait:
        #write.table(fits.criteria.raw.master, paste("Raw.P_",trait,"_Genes.pl.min.", search.range,"bp.of.GWAS.SNPs_QTL",cSI,".txt", sep = ""), sep = "\t", quote = FALSE)
        #write.table(fits.criteria.FDR.master, paste("FDR.P_",trait,"_Genes.pl.min.", search.range,"bp.of.GWAS.SNPs_QTL",cSI,".txt", sep = ""), sep = "\t", quote = FALSE)              
      
        }#end if("search.range")
      
      if(pvals.to.output == "1gene"){ # end search.range loop, start 1-gene
        print(paste("Now processing, for 1gene, trait ",trait,sep=''))
        
        #Read in proximal gene lists produced by "Identify_genes_proximal_to_RMIP_SNPs.r"
        proximal.gene.data = read.table(paste(home.dir,proximal.genes.dir,"Genes.proximal.to.RMIP.SNPs.for.",trait,".txt", sep=""), sep ='\t',head=TRUE,stringsAsFactors = FALSE)

        #Iterate through RMIP SNPs for this trait-SI pair
        for(SNP.index in (1:nrow(RMIP.SNPs.this.trait_cSI))){
          #Determine RMIP SNP chr and pos from RMIP common S.I. key subset, columns 1 and 2
          this.SNP.chr = RMIP.SNPs.this.trait_cSI[SNP.index,chr.col.common.key]
          this.SNP.pos = RMIP.SNPs.this.trait_cSI[SNP.index,pos.col.common.key]
          
          #Extract proximal genes corresponding to this RMIP SNP
          proximal.genes.to.this.SNP = proximal.gene.data[which(
            proximal.gene.data[,2]==this.SNP.chr & 
            proximal.gene.data[,3]==this.SNP.pos),]
          #Remove NA values (indicate no exact hit)
          proximal.genes.to.this.SNP = proximal.genes.to.this.SNP[which(!(is.na(proximal.genes.to.this.SNP[,6]))),]
          
          #Extract raw p-values for these proximal genes and append them as column to proximal gene list for this RMIP SNP
          proximal.summary.this.SNP = merge(proximal.genes.to.this.SNP,the.raw.P.values[,c(1:9)],by.x = 6, by.y = 1)
          
          #Append proximal genes with p values to master file for this common S.I.
          proximal.summary.master = rbind(proximal.summary.master,proximal.summary.this.SNP)
        } #end loop through RMIP SNPs

        #####If want to output by single trait:
        #write.table(proximal.summary.master,paste("Raw.P.values_",trait,"_Genes.hit.or.adjac.to.GWAS.SNPs_QTL",cSI,".txt",sep=''),sep='\t',row.names=FALSE)   
      } #end if("1gene")
    }#end loop indicating non-zero number of SNPs to process
} #End loop through traits

  #Output master files for this common S.I.
  if(pvals.to.output == 'all.in.interval'){
    #Master data frame is already sorted by trait, then starting (left) bp position; no need to sort.
    write.table(FDR.pval.this.cSI_allTraits,paste("FDR.P.values_all.in.interval_QTL",cSI,".txt",sep=''),sep='\t',row.names=FALSE)
  }# end all.in.interval
  if(pvals.to.output == "search.range"){
    if(!(is.null(fits.criteria.raw.master))){
    raw.search.range.sorted = fits.criteria.raw.master[with(fits.criteria.raw.master,order(fits.criteria.raw.master$trait,fits.criteria.raw.master$Start_bp)),]
    write.table(fits.criteria.raw.master, paste("Raw.P_Genes.", search.range,"bp.of.GWAS.SNPs_QTL",cSI,".txt", sep = ""), sep = "\t", quote = FALSE,row.names=FALSE)
    #FDR.search.range.sorted = fits.criteria.FDR.master[order("trait","Chr","Start_bp"),]
    #write.table(fits.criteria.FDR.master, paste("FDR.P_Genes.", search.range,"bp.of.GWAS.SNPs_QTL",cSI,".txt", sep = ""), sep = "\t", quote = FALSE)  
    
    #Export the FDR-adjusted P-values, adjusted separately for each gene for the search space of its respective RMIP SNP
    write.table(FDR.adjusted.P.values.all.SNPs.this.QTL, paste(home.dir, fdr.adjusted.for.100kb.dir,"FDR.Adjusted.P.values_by_100kb_QTL",cSI,".txt", sep = ""), 
                quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
    }#end if is not null
    }# end search.range
  if(pvals.to.output == "1gene"){
    if(!(is.null(proximal.summary.master))){
      proximal.summary.master.sorted <- proximal.summary.master[with(proximal.summary.master,order(proximal.summary.master$trait,proximal.summary.master$SNP.position)),]
      write.table(proximal.summary.master.sorted,paste("Raw.P.values_Genes.hit.or.adjac.to.GWAS.SNPs_QTL",cSI,".txt",sep=''),sep='\t',row.names=FALSE)
    }#end if is not null
  }# end 1gene
} #End loop through common S.I.s
