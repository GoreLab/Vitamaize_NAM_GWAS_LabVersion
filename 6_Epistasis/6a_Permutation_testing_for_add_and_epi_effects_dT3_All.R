rm(list = ls())
#trait residuals corrected for Y1 or major gene
#code evaluates two models
#strategy 1 generates tested model with all additive effects and possible epistatic effects, generates permutation threshold using lm, selects terms according to threshold, refits final model
#strategy 2 generates tested model with all additive effects and possible epistatic effects, tests with stepwise regression using AIC as model selection criterion

###Required files:
#(1) GBS genotype files
#(2) Transformed BLUE files
#(3) BLUE residuals from Y1 model
#(4) JL model results

#######################################################################################################################################################
#USER SPECIFIED INPUTS
setwd("/Users/anybody/Documents/Toco_NAM/")
home.dir <- getwd()
geno.dir = paste(home.dir,"/Geno.Pheno_Inputs/",sep='')
pheno.dir = geno.dir
dT3.new.path = paste(home.dir,"/Methods/dT3_removeExtremeVal_test/",sep='')
dir.for.modified.JL.results <- paste(home.dir,"/JL/Toco_FinalModels_postMCCorrPostRescan/",sep='')
pairwise.epistasis.results.dir <- pheno.dir
permut.dir <- paste(home.dir,"/Epistasis/",sep='')
final.combined.model.fit.and.PVE.dir <- permut.dir

PVE.byFamily = TRUE

library(multtest)

#Read in the genotypic data
genotypes <- read.table(paste(geno.dir,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#select traits, if there are any file name changes specific to MC traits, list MC trait vector here, if not leave trait.collinear blank
trait = c("gT","gT3","aT","aT3","dT","dT3","PC8","totalT","totalT3","totalTocochrs")
trait.collinear <- c("aT3","dT3","gT3","totalT3","totalTocochrs")   

#Specify types of analyses to be run; set those desired to TRUE
generate.model = TRUE                          #Generate complete model with all additive and possible epistatic terms except those including major gene
generate.permutations = TRUE                 #Requires generate.model to be true; randomize residuals by family across NAM; generate designated number of permutations; identify thresholds with desired typeI error rate
                                                ###NOTE REmove comment from bracket in line 206 to generate permutations and loop around traits
test.epi.interactions.in.final.model = TRUE   #Requires generate.model to be true; Obtain pval for all interaction terms when tested singly in additive model (minus major gene); apply threshold and calculate PVE of final model terms
test.stepwise.regression.model = FALSE         #Requires generate.model to be true; identify optimal model from all terms (minus major gene) using AIC model selction criterion

#Specify number of permutations to conduct in generate.permutations
number.of.permuted.phenos <- 1000
permutation.testing.results <- NULL

#Specify typeI error rate at which to conduct final model fit
select.type1.error.rate <- "0.01"           #choices can be however many user allows in the generate.permutations section; this variable will be called in to generate final model fit in test.epi.interactions.in.final.model

count_AltAllele = function(colIndex){
  m=colIndex
  print(m)
  minor.counts.vector = as.vector(tapply(genoScores[,m], genoScores[,2],function(scoreSet){
    s=as.numeric(scoreSet)
    #return(1*length(which(s<=0.5 & s<1.5)+2*length(which(s>=1.5 & s<=2)))) #if trying to use distance imputed values--leads to MAFs of 0.2 to 0.3
    return(1*length(which(s==1))+2*length(which(s==2)))
  }))
  altAF = minor.counts.vector/(2*nf.vector)
  print(altAF)
  #print(maf)
  pooled_altAF = (1/N)*crossprod(nf.vector, altAF)
  print(paste("Pooled Frequency of alternate allele ignoring distance-imputed values is",pooled_altAF))
}
######################################################################################################################################################
#For loop through the traits
for(i in trait){

if(generate.model == TRUE){

  #Read in the trait under study, probably don't need this data set
  if(i == "dT3"){
    pheno.data = read.table(paste(dT3.new.path,"Box-Cox_TRANS_filesforjlremoveextremeval_dt3/BLUEs_No_Outliers_Transformed_dT3.only.Redone.txt",sep=''),head=TRUE)
    pheno.data$dT3 = as.numeric(pheno.data$dT3)
  }else{
  pheno.data = read.table(paste(pheno.dir,"BLUEs_No_Outliers_Transformed_all_for_TASSEL_",i,".txt", sep = ""), head = TRUE)
  }
  
  colnames(pheno.data) <- c("Geno_Code", i)

  #Read in the results from JL analysis
    if(i=="dT3"){
      TASSEL.model.results <- read.table(paste(dT3.new.path,"new.trans_new.perm_FINAL/JL_fromTrans_dT3_Redone_MC.corrected_Tassel3fromSF_forR.txt",sep=''), head = TRUE)
    }else if(i %in% trait.collinear){
    TASSEL.model.results <- read.table(paste(dir.for.modified.JL.results,"MC_corrected_",i,"_postRescan_R.formatted.txt"  , sep = "") , head = TRUE)
    }else{
    TASSEL.model.results <- read.table(paste(dir.for.modified.JL.results,i,"_model_3fromSF_0.01_Final_R.formatted.txt", sep = "") , head = TRUE)
    }
  
  #Parse out all of the SNPs that were selected in the model
  geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,4]),-c(2:5)]
  chr.reduced <- TASSEL.model.results[,2][-1]
  cM.pos.reduced <- TASSEL.model.results[,3][-1]
  geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
  colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))

  #Get the appropriate SNPs and merge the phenotypic and genotypic data together.
  geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

  #Add a population column
  geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
  colnames(geno.and.pheno)[2] <- "pop"

  #Get the model order of the model in R to be the same as that in TASSEL.
  seq <- (seq(1:nrow(TASSEL.model.results))-1)[-1]
  model.order <- cbind(seq, TASSEL.model.results[-1,c(2,3,4)])
  #Sort by chromosome and bp
  model.order <- model.order[order(model.order[,3]),]
  model.order <- model.order[order(model.order[,2]),]
  model.order <- cbind(model.order, seq(1:nrow(model.order)))
  #Sort by the first column of marker order so that the markers will be put into the model in the same order as they were selected.
  model.order <- model.order[order(model.order[,1]),]
  index <- model.order[,5]+3 #3 is added so because the first SNP is on the fourth column

  #generate additive model "base.model"
  base.model <- paste(i, " ~ pop",sep = "")
  for(k in index){
    base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")
  } #end for(k in 4:ncol(geno.and.resid))

  #obtain vector of additive effects
  the.markers <- as.character(TASSEL.model.results[-1,4])

  #generate vector of all pairwise combinations, for j additive effects there should be !(j-1) interaction terms
    all.interaction.terms.vector <- NULL
   for(j in 1:(length(the.markers)-1)){
    for(k in (j+1):length(the.markers)){
        single.interaction.term <- paste("pop:",the.markers[j],":",the.markers[k], sep = "")
        all.interaction.terms.vector <- c(all.interaction.terms.vector, single.interaction.term)
     } #end k loop
    } #end j loop

  #format all interaction vector to fit into model
  the.interaction.term.model.component <- NULL
  combine.model.components <- NULL
  for (m in all.interaction.terms.vector[1:(length(all.interaction.terms.vector)-1)]){
    the.interaction.term.model.component <- paste(m,"+")
    combine.model.components <- paste(combine.model.components, the.interaction.term.model.component)
    }  # end adding interaction terms with the exception of the last term in the vector

    #add last term in vector
    final.interaction.term.model.component <- paste(combine.model.components, all.interaction.terms.vector[length(all.interaction.terms.vector)])

  #Generate complete model with all additive and possible epistatic terms
      #the.tested.model <- paste(base.model, "+",final.interaction.term.model.component ,sep = "")
      the.tested.model <- base.model
   } #end generate model

 #######################################################################################################################################
  if(generate.permutations == TRUE){
  
  type1.error.rate.01 <- 0.01
  type1.error.rate.05 <- 0.05
  
 #Generate permuted residuals within populations
 #Test each set with the.tested.model
 #Compile p values from each set

 all.permutation.results <- NULL
 for (p in 1:number.of.permuted.phenos){

 complete.this.NAM.permutation <- NULL
  for (u in unique(geno.and.pheno[,2])){    #permute by pop
  ##########Add family permutations
    pop.select <- geno.and.pheno[which(geno.and.pheno[,2] == u),]
    pop.select.pheno <- pop.select[,3]
    pop.shuffle <- sample(pop.select.pheno, length(pop.select.pheno))
    permute.pheno.data <- cbind(as.data.frame(pop.select[,1:2]), as.data.frame(pop.shuffle),as.data.frame(pop.select[,4:length(pop.select)]))
    complete.this.NAM.permutation <- rbind (complete.this.NAM.permutation, permute.pheno.data)
    } # end u pop permutations
    colnames(complete.this.NAM.permutation)[3] <- i

          #Fit the model with all terms
          JL.permutation.model <- lm(paste(the.tested.model, sep = ""), data = complete.this.NAM.permutation)
          #ANOVA.two.way.interaction.this.model <-  anova(JL.model)

          ANOVA.two.way.interaction.this.model <-  as.data.frame(anova(JL.permutation.model)[2:(nrow( anova(JL.permutation.model))-1),])

          permutation.vector  <- rep(p, nrow(ANOVA.two.way.interaction.this.model))
          results.this.model <- as.data.frame(cbind(rownames(ANOVA.two.way.interaction.this.model), ANOVA.two.way.interaction.this.model, permutation.vector))
          all.permutation.results <- rbind(all.permutation.results, results.this.model)
 } #end p residuals

 #Combine all p values and sort
  all.permuted.pvalues <- sort(all.permutation.results[,6])
  
 #Identify thresholds
  threshold.01 <- all.permuted.pvalues[round(length(all.permuted.pvalues)*type1.error.rate.01)]
  threshold.05 <- all.permuted.pvalues[round(length(all.permuted.pvalues)*type1.error.rate.05)]
  
 #Append data for all traits
  print(paste("Generating permutation thresholds for trait ", i, " using ",  number.of.permuted.phenos, " permutations", sep = ""))
  data.keep.01 <- c(i, type1.error.rate.01, number.of.permuted.phenos, threshold.01)
  data.keep.05 <- c(i, type1.error.rate.05, number.of.permuted.phenos, threshold.05)
  permutation.testing.results <- rbind(permutation.testing.results, data.keep.01, data.keep.05)
   #} #end i traits CHD commented out; seems to be ending analysis by trait too early.
  write.table(all.permuted.pvalues, paste(permut.dir, "Permuted_data_from_tested_additive_model_for_",i,".txt", sep=""), sep = "\t", row.names = FALSE)
  colnames(permutation.testing.results) <- c("Trait", "TypeI_Error_Rate", "No.Permutations", "Threshold")
  write.table(permutation.testing.results, paste(permut.dir,"Permutation_thresholds_for_single_epistatic_interaction_testing_",i,".txt", sep = ""), sep = "\t", row.names = FALSE)
  } #end generate.permutations

  #######################################################################################################################
  if(test.epi.interactions.in.final.model == TRUE){
  
  #requires base model from generate.model
  #do not use the.tested.model from above, as it contains all additive and possible epistatic terms
  #generate models with single epistatic terms
  
  #For loop through all possible pairwise combinations of models
  for(j in 1:(length(the.markers)-1)){
    for(k in (j+1):length(the.markers)){
      the.interaction.term <- paste("pop:",the.markers[j],":",the.markers[k], sep = "")
      the.tested.single.model <- paste(base.model, "+",the.interaction.term ,sep = "")

      #Fit the base model and the j,kth interaction term
      JL.model <- lm(paste(the.tested.single.model, sep = ""), data = geno.and.pheno)
      ANOVA.two.way.interaction.this.model <-  anova(JL.model)[nrow( anova(JL.model))-1,]

       #######Some code may heve been omitted at this point
       #######Check original version if this generates an error

      #Append the corresponding row of the ANOVA table to the results
      if((j==1)&(k==2)){
        results <- ANOVA.two.way.interaction.this.model
      }else{
        results <- rbind(results, ANOVA.two.way.interaction.this.model)
      }#end if((j==1)&(k==2))

    } #end for(k in (j+1):length(the.markers))
  }#End for(j in 1:(length(the.markers)-1))

  #Sort results smallest to largest P-values and add FDR adjustment
  results <- cbind(rownames(results), results)
  results <- results[order(results[,ncol(results)]),]
  FDR.adj.pval <- p.adjust(results[,6], method = "BH")
  results <- cbind(results, FDR.adj.pval)
  colnames(results)[1] <- "Term(NOTE:A_separate_model_was_fitted_for_each_two-way_interaction_term)"
  colnames(results)[7] <- "FDR.adj.pval"
  
  #write.table(results, paste(home.dir, pairwise.epistasis.results.dir,"Two.way.epistasis.models.for.",i,".NAM.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
 
    #identify interaction terms from single term testing above that pass the specified threshold and refit final model
    #call in threshold table, specify by trait i and error rate (above)
    ##recorded.thresholds <- read.table(paste(home.dir, permut.dir, "Permutation_thresholds_for_single_epistatic_interaction_testing_combined_alpha.txt", sep = ""), head = TRUE)
  recorded.thresholds <- read.table(paste(permut.dir, "Permutation_thresholds_for_single_epistatic_interaction_testing_",i,".txt", sep = ""), head = TRUE)
    select.threshold <- recorded.thresholds[which((recorded.thresholds[,1] == i) & (recorded.thresholds[,2] == select.type1.error.rate)),4]
    
    #identify interaction terms that are significant (pval, not FDR adj Pval) at or below the threshold
    select.interactions <- as.vector(results[which(results[,6] <= select.threshold), 1])

    #add these terms to base model
    if(length(select.interactions) == "0"){
    the.interaction.terms <- NULL}else{
      the.interaction.terms <- NULL
      for(p in 1:length(select.interactions)) {
      the.interaction.terms <- paste(the.interaction.terms,"+",select.interactions[p],sep = "")
      } #end p terms added to model
    } #end if loop accounting for no interaction terms
    the.tested.final.model <- paste(base.model,the.interaction.terms ,sep = "")
    
    #Fit the all inclusive model
    JL.model <- lm(paste(the.tested.final.model, sep = ""), data = geno.and.pheno)
    ANOVA.entire.model <- anova(JL.model)
    
    #generate PVE
    ########Generating variance and sample size vectors
    ####Extract the trait value, RIL ID, and family
    pheno.data.temp <- geno.and.pheno[,1:3]

    #Remove all missing phenotypic data
    if(length(which(is.na(pheno.data.temp[,3])))>0) pheno.data.temp <- pheno.data.temp[-which(is.na(pheno.data.temp[,3])),]
    print("1")
    ####Calculate the following summary statistics:
    #Variance within each family
    Vf.vector <- as.vector(tapply(pheno.data.temp[,3], pheno.data.temp[,2], var))
    print("2")
    #Sample size within each family
    nf.vector <- as.vector(tapply(pheno.data.temp[,3], pheno.data.temp[,2], length))

        #pop matrix used to exclude pops with incomplete complement of interactions in epistatic marker effects
        pop.vector <- as.vector(unique(pheno.data.temp[,2]))
        pop.nf.matrix <- cbind(pop.vector, nf.vector)

    #Overall sample size
    N <- sum(nf.vector)
    nf.vector.by.N <- nf.vector/N

    den <- (1/N)*crossprod(nf.vector, Vf.vector)
    print("3")
  
      #############generate effect estimates for all terms in model
  all.terms <- as.vector((rownames(ANOVA.entire.model)))
  all.terms <- all.terms[2:(length(all.terms)-1)]
  all.terms.no.pop <- substr(all.terms, start = 5, stop = 1000)

  term.vector <- NULL
  PVE.vector <- NULL
  
  pop.Effect_allInteraxns = NULL
  
    #Obtain marker effect estimates for additive terms
     for(r in all.terms.no.pop){
    #epistatic.marker.effects <- summary(JL.model)$coefficients[which(substr(rownames(summary(JL.model)$coefficients), start = 7, stop = 1000) == substr(the.interaction.term,start = 5, stop = 1000)) ,]
     marker.effects <- summary(JL.model)$coefficients[which(substr(rownames(summary(JL.model)$coefficients), start = 7, stop = 1000) == r) ,]
     thisTerm = strsplit(r,split=":")
     if(length(thisTerm[[1]]) > 1){
        thisTermEpistatic = TRUE
        marker1= thisTerm[[1]][1]
        marker2 = thisTerm[[1]][2]
     }else{thisTermEpistatic = FALSE}
     
       #obtain pops represented in epistatic marker effects - CBK ADDED
        pop.ID <- substr(rownames(marker.effects), start=4, stop=5)

       #reduce pop.nf.matrix to only those pops represented in epistatic effects - CBK ADDED
           pop.match <- pop.nf.matrix[which(pop.nf.matrix[,1] %in% pop.ID),]

       #obtain new nf.vector with appropriate number of pops - CBK ADDED
        nf.vector.mod <- as.numeric(as.vector(pop.match[,2]))
        nf.vector <- as.vector(nf.vector.mod)

    if (PVE.byFamily == FALSE){
      sq.term <- marker.effects[,1]^2
      num.left.term  <- ((2/N)*(crossprod(nf.vector,sq.term)))
      num.right.term <- ((1/N)*(crossprod(nf.vector,marker.effects[,1])))^2
      num <- num.left.term - num.right.term
      PVE.by.term <- num/den
    }else if (PVE.byFamily == TRUE){
      colnames(geno.and.pheno)[1] = "Taxa"
      if(thisTermEpistatic == TRUE){
        thisGeno.with.Taxa = geno.and.pheno[c("Taxa",marker1,marker2)]
      }
      if(thisTermEpistatic == FALSE){
        thisGeno.with.Taxa = geno.and.pheno[c("Taxa",r)]
      }
      genoScores = merge(pheno.data.temp,thisGeno.with.Taxa,by.x=1,by.y = 1) #append the genos just for this marker to line,pop,pheno cols           
      genoScores = subset(genoScores,genoScores[,2] %in% pop.match[,1])
      
      #MAF for this SNP (or 1st SNP in interaction term)
      count_AltAllele(4)
      #MAF for 2nd SNP in interaction term (when applicable)
      if(thisTermEpistatic==TRUE){
        count_AltAllele(5)
      }
      
      #Beginning of actual PVE calculation
      if(thisTermEpistatic==FALSE){
        genoScores[,4] = genoScores[,4] - 1 #need to have -1,0,1 coding so that homozyg minor allele has effect of -n
        genoScores.by.Effects.vector = genoScores[,4] 
      }
      
      if(thisTermEpistatic==TRUE){
        genoScores[,c(4,5)] = genoScores[,c(4,5)] - 1 #need to have -1,0,1 coding so that homozyg minor allele has effect of -n
      
        #Frequencies of each genotype class for the interaction
        classFreqs = ftable(genoScores[,c(4,5)],row.vars=marker1,col.vars=marker2)
        classFreqs_df = as.data.frame(classFreqs)
        order.classes = order(classFreqs_df$Freq, decreasing=T)
        View(classFreqs_df[order.classes,])
      
        compositeGeno = genoScores[,4]*genoScores[,5]
        genoScores.by.Effects.vector = compositeGeno
      }
      
      pop.Effect = cbind(substr(rownames(marker.effects), start=4, stop=5),marker.effects[,1])
      pop.Effect_allInteraxns = rbind(pop.Effect_allInteraxns,pop.Effect)
      
      for (this.pop in pop.match[,1]){
        #print(this.pop)
        Effect = as.numeric(pop.Effect[which(pop.Effect[,1]==this.pop),2])
        for (this.row in 1:nrow(genoScores)){
          #print(paste("Value of effect*genoScore vector before multiplying (should just be geno score):",genoScores.by.Effects.vector[this.row],sep=''))
          
          #Calculate effect * marker score, *-1 because effect estimate  (including sign) is for B73 reference allele, which is now coded as -1
          if(genoScores[this.row,2]==this.pop){genoScores.by.Effects.vector[this.row]=-1*Effect*genoScores.by.Effects.vector[this.row]}
          
          #print(paste("Value of effect*genoScore vector after multiplying:",genoScores.by.Effects.vector[this.row],sep=""))
        }
      }
      pop_Scores.by.Effects = cbind(genoScores[,2],genoScores.by.Effects.vector)
      Vg.vector = tapply(pop_Scores.by.Effects[,2],pop_Scores.by.Effects[,1],var)
      Vg.den <- (1/N)*crossprod(nf.vector, Vg.vector) 
      PVE.by.term <- Vg.den/den
    } 

    pop.term <- paste("pop:", r, sep="")
    term.vector <- c(term.vector, pop.term)
    PVE.vector <- c(PVE.vector, PVE.by.term)

     } #end r loop for interaction term PVE calculation

  print("4")
  
  write.table(pop.Effect_allInteraxns,paste(final.combined.model.fit.and.PVE.dir,"/Pop.by.Marker_Effect.Estims_Epist_FinalJLModel_",i,".txt",sep=''),sep="\t")
  
  PVE.matrix <- cbind(term.vector, PVE.vector)
  output <- merge(ANOVA.entire.model, PVE.matrix, by.x = "row.names", by.y = "term.vector")

   write.table(output, paste(final.combined.model.fit.and.PVE.dir,"Add.and_epi.JL.model.plus.PVE.at.alpha.",select.type1.error.rate,"for.",i,".NAM_chd_PVEbyFam.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
  
  #} #end i traits
  } # end test.epi.interactions.in.final.model

  ######################################################################################################################
  if(test.stepwise.regression.model == TRUE){
  
  #test model with stepwise regression
  JL.model <- lm(paste(the.tested.model, sep = ""), data = geno.and.pheno)
  step.reg.backward <- step(JL.model, direction = "backward")
  #step.reg.forward <- step(JL.model, direction = "forward")
  #step.reg.both <- step(JL.model, direction = "both")

  #Obtain ANOVA model fit and PVE from optimal model from stepwise selection
    #generate new tested model equivalent with lowest AIC
    #add PVE code
    
    } #end test.stepwise.regression.model
    
} #end trait i