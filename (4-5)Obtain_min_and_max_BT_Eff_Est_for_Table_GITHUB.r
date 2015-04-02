rm(list = ls())

###Required files:
### (1) Modified tabular summary file from script (4-3)
### (2) BACKTRANSFORMED marker effect estimates from script (2-3)



the.extractinator.for.Supp.Tab.2 <- function(){
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
  founder.names <- as.data.frame(c("B97", "CML103", "CML228", "CML247", "CML277", "CML322", "CML333", "CML52", 
                                   "CML69", "HP301", "Il14H", "KI11", "KI3", "KY21", "M162W", "M37W", "MO18W", 
                                   "MS71", "NC350", "NC358", "OH43", "OH7B", "P39", "TX303", "TZI8"))
  NAM.pops <- cbind(pop.seq, founder.names)
  colnames(NAM.pops) <- c("Pop.num", "Pop.Founders")
  
  #Obtain the populationxmarker effect estimates for this QTL
  
  effect.estimates.from.entire.trait <- read.table(paste(home.dir,location.of.effect.estimates,
                                         "NonStd.Pop.by.Marker.Effect.Est.from.R.B73centered.SI01.",trait,".txt", sep = ""),head = TRUE)
  
  
 # multicollinear.traits <- c("BCAR_RUV", "TOTCAR_RUV")
  #if (trait %in% multicollinear.traits){
  #effect.estimates.from.entire.trait <- read.table(paste(home.dir,location.of.effect.estimates,
  #                                       "NonStd.Pop.by.Marker.Effect.Est.from.R.B73centered.",trait,".no.col.txt", sep = ""),head = TRUE)
  #}else{
  #effect.estimates.from.entire.trait <- read.table(paste(home.dir,location.of.effect.estimates,
  #                                       "NonStd.Pop.by.Marker.Effect.Est.from.R.B73centered.",trait,".txt", sep = ""),head = TRUE)
  #}
  
  markerID <-  substr(effect.estimates.from.entire.trait[,1], 7,10000)
  popID <- substr(effect.estimates.from.entire.trait[,1], 1,5)
  
  
  bt.effect.estimates.from.QTL <-data.frame(as.character(popID[which(markerID == QTL)]),
                                  as.numeric(effect.estimates.from.entire.trait[which(markerID == QTL), ncol(effect.estimates.from.entire.trait)])
                                  )
  colnames(bt.effect.estimates.from.QTL) <- c("Population","Back.Trans.Est")
  
  bt.effect.estimates.from.QTL <- merge(bt.effect.estimates.from.QTL, NAM.pops, by.x = "Population", by.y = "Pop.num")
  

  return(bt.effect.estimates.from.QTL)
}#end the.extractinator.for.Supp.Tab.2 



##########################################################################
##########################################################################
##########################################################################







#Set up all of the paths and the input and output file names
setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
home.dir <- getwd()
location.of.effect.estimates <- "\\(12)Standardizing Allelic Effects\\"
location.of.GWAS.results <- "\\(10)GWAS Analysis\\RUV and MSS folders - original 25fam analyses\\"
output.dir <- "\\(16)Generating Robust Files for Group Review\\(16c)Pipeline for screening candidate RMIP\\Combining_JL_GWAS_FPKM_Workspace\\"
#ratios.or.no.ratios <- "No_Ratios"  #Options are "No_Ratios" or "with_Ratios"
input.file.name <- "Tab_Sum_Carot_with_Common_SI_Info_updated_20140617_left_bound.txt"


#Read in the tabular summary
setwd(paste(home.dir,"\\(16)Generating Robust Files for Group Review\\Generating overlapping support intervals\\",sep = ""))
tabular.summary <- read.table(paste(input.file.name,sep = "") ,head = TRUE)

minimum.effect.est <- NULL
minimum.effect.est.family <- NULL
maximum.effect.est <- NULL
maximum.effect.est.family <- NULL

for(i in 1:nrow(tabular.summary)){#For loop through the rows of the summary
  #Run the.extractinator, which will output the back-transformed alleleic effect estimates
  trait <- as.character(tabular.summary[i,1]) #"DT"
  QTL <- as.character(tabular.summary[i,6]) #"S_204809285"
  Gene <- NA
  Gene.GRZM.ID <- NA  
  
  
  bt.effect.estimates.from.QTL <- the.extractinator.for.Supp.Tab.2()
  
  #Record the minimum allelic effect estimate into one vector, and the corresponding family
  minimum.effect.est <- c(minimum.effect.est, bt.effect.estimates.from.QTL[which(bt.effect.estimates.from.QTL[,2] == 
                          min(bt.effect.estimates.from.QTL[,2])), 2])

  minimum.effect.est.family <- c(minimum.effect.est.family, paste("B73x",bt.effect.estimates.from.QTL[which(bt.effect.estimates.from.QTL[,2] == 
                          min(bt.effect.estimates.from.QTL[,2])), 3],sep = ""))  

  #Record the maximum allelic effect estimate into one vector, and the corresponding family
  maximum.effect.est <- c(maximum.effect.est, bt.effect.estimates.from.QTL[which(bt.effect.estimates.from.QTL[,2] == 
                          max(bt.effect.estimates.from.QTL[,2])), 2])

  maximum.effect.est.family <- c(maximum.effect.est.family, paste("B73x",bt.effect.estimates.from.QTL[which(bt.effect.estimates.from.QTL[,2] == 
                          max(bt.effect.estimates.from.QTL[,2])), 3],sep = ""))  
}#End for loop through the rows of the summary

#Append the minimum and maximum allelic effect estimate value vectors, and the corresponding family

tabular.summary <- cbind(tabular.summary, minimum.effect.est, minimum.effect.est.family, 
                         maximum.effect.est, maximum.effect.est.family)


#Export the resutls
write.table(tabular.summary,  paste("Tabular.summary.for.Supp.Tab.2.Carots.bt.eff.ests.txt", sep = ""),
                                    quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)


