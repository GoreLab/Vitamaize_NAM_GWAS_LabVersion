require("ggplot2")
require("ggthemes")
require("Cairo")
require("grid")
require("gridExtra")
require("scales")
require("zoo")
require("reshape2")
require("ggrepel")


##################################
##################################
#
# LOCAL FUNCTIONS
#
##################################
##################################


theme_blank = function(base_size = 12, base_family = "Helvetica"){
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
        rect = element_blank(),
        line = element_blank(),
        text = element_blank(),
        axis.text = element_text(margin = unit(0, "lines"))
        )
}

plot_blank = function(){

    empty = ggplot() +
        theme_blank() +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        theme(axis.text = element_blank(), panel.margin = unit(0, "mm")) +
        labs(title=NULL, x=NULL,y=NULL)

    return(empty)
}


rescaleMb = function(bp){
    scale = 1000000
    l = round(bp / scale, digits = 2)
    comma(l)
}






##################################
##################################
#
# PARAMETERS SETUP
#
##################################
##################################


##########################
# Read user parameters from environment
# and overwrite when/if necessary

# read relevant environment variables
env.var = as.list(
    Sys.getenv(
        x = c(
            "fig_d",
            "output_f",
            "output_file_type",
            "ld_f",
            "expr_f",
            "expr_sig_f",
            "rmip_f",
            "genes_f",
            "key_f",
            "search_space",
            "bin_size",
            "pval_cutoff",
            "ribbon_type",
            "label_type"
            ),
        unset = NA,
        names = TRUE
        )
    );

# store I/O parameters
out.dir = as.character(env.var$fig_d)
out.file = as.character(env.var$output_f)
out.type = as.character(env.var$output_file_type)
ld.file = as.character(env.var$ld_f)
expr.file = as.character(env.var$expr_f)
expr.sig.file = as.character(env.var$expr_sig_f)
gwas.file = as.character(env.var$rmip_f)
genes.file = as.character(env.var$genes_f)
key.file = as.character(env.var$key_f)


# store data transform parameters
kSearchSpace = as.integer(env.var$search_space)
kBinSize = as.integer(env.var$bin_size)
kExprPvalCutoff = as.numeric(env.var$pval_cutoff)

# store visualization parameters

ribbon.type = as.character(env.var$ribbon_type)
label.type = as.character(env.var$label_type)

############################
#Modify environment vars as needed

#out.dir = as.character(env.var$fig_d)
#out.file = as.character(env.var$output_f)
#out.type = as.character(env.var$output_file_type)
#ld.file = as.character(env.var$transformed_ld_file)
#expr.file = as.character(env.var$transformed_expr_file)
#expr.sig.file = as.character(env.var$transformed_expr_sig_file)
#gwas.file = as.character(env.var$transformed_gwas_file)
#genes.file = as.character(env.var$transformed_genes_file)
#key.file = as.character(env.var$transformed_key_file)
#out.dir = as.character(env.var$output_dir)

#kExprPvalCutoff = 0.05
#kSearchSpace = 100000
#kBinSize = 2000
############################


###########################
# Other params
#--------------------------
kMinRsq = 1
kDensityAdjust = 0.01
kGeneArrowProp = 0.2
###########################


###########################
# Variant params
#--------------------------
ribbon.type.default = 'flat'
label.type.default = 'arrow'

# uncomment to change/overwrite

ribbon.type = 'flat'
#ribbon.type = 'hist'

label.type = 'arrow'
#label.type = 'text'
#label.type = 'off'

# defaults for variants
if(! exists("ribbon.type")){
    ribbon.type = ribbon.type.default
}

if(! exists("label.type")){
    label.type = label.type.default
}



###########################



##########################
# Graphic parameters
kBorderColour = "black"
kBorderSize = 0.3

kTitleSize = 10
kTitleColour = "black"
kTitleFace = "bold"

kSubTitleSize = 8
kSubTitleColour = "black"
kSubTitleFace = "bold"

kAxisLabelSize = 6
kAxisLabelColour = "black"
kAxisLabelFace = "plain"

kAxisTextSize = 6
kAxisTextColour = "black"
kAxisTextFace = "plain"
##########################


##########################
# Ordering, labeling, colours and shapes
##########################

#Set sample order for expression correlation
samples.order = c(
    "12DAP",
    "16DAP",
    "20DAP",
    "24DAP",
    "30DAP",
    "36DAP"
    )

# Set traits order, labels and other formatting parameters
trait.order = c(
    "aT",
    "dT",
    "gT",
    "aT3",
    "dT3",
    "gT3",
    "totalT",
    "totalT3",
    "totalTocochrs",
    "PC8"
    )

trait.group.order = c(
    "T",
    "T3",
    "Sum",
    "PC"
    )


#Colours based on compound class (bond saturation):
# If trait ends with T (tocopherol): "#fec44f" 
# If trait ends with T3 (tocotrienol): "#d95f0e"
# If total T: "#6baed6"
# If total T3: "#3182bd"
# If total T+T3: "#08519c"
# If PC8: "#762a83"

#Shapes based on isoforms (placement of methyl groups):
#  If trait starts with 'a' (representing alpha isoform): square
#  If 'd' (delta): rhomboid
#  If 'g' (gamma): circle
#  If 't' (total--these are derived, sum traits): upwards triangle
#  If PC8: downwards triangle


trait.format = data.frame(
    Trait = trait.order,
    TraitGroup = c(
        "T",
        "T",
        "T",
        "T3",
        "T3",
        "T3",
        "Sum",
        "Sum",
        "Sum",
        "PC"
        ),
    TraitLabel = c(
        "alpha*T",
        "delta*T",
        "gamma*T",
        "alpha*T*scriptstyle(3)",
        "delta*T*scriptstyle(3)",
        "gamma*T*scriptstyle(3)",
        "Sigma*(T)",
        "Sigma*(T*scriptstyle(3))",
        "Sigma*(list(T, T*scriptstyle(3)))",
        "PC*8"
        ),
    TraitColour = c(
        "#fec44f",
        "#fec44f",
        "#fec44f",
        "#d95f0e",
        "#d95f0e",
        "#d95f0e",
        "#6baed6",
        "#3182bd",
        "#08519c",
        "#762a83"        
        ),
    TraitShape = c(
        22,
        23,
        21,
        22,
        23,
        21,
        24,
        24,
        24,
        25
        )
)
##########################










##################################
##################################
#
# DATA PREP
#
##################################
##################################


##########################
# read data from file
ld.SNPs = read.table(ld.file, header = TRUE, sep = "\t",strip.white = TRUE, comment.char = "", quote = "", na.strings = c("?", ""));
GWAS.hits = read.table(gwas.file, header = TRUE, sep = "\t",strip.white = TRUE, comment.char = "", quote = "", na.strings = c("?", ""));
expr.corr = read.table(expr.file, header = TRUE, sep = "\t",strip.white = TRUE, comment.char = "", quote = "", na.strings = c("?", ""));
expr.corr.pval = read.table(expr.sig.file, header = TRUE, sep = "\t",strip.white = TRUE, comment.char = "", quote = "", na.strings = c("?", ""));
genes.coord = read.table(genes.file, header = TRUE, sep = "\t",strip.white = TRUE, comment.char = "", quote = "", na.strings = c("?", ""));
linking.key = read.table(key.file, header = TRUE, sep = "\t",strip.white = TRUE, comment.char = "", quote = "", na.strings = c("?", ""));
##########################





# Add relevant key columns to input data for linking/selecting across datasets
linking.key$QTL_No = factor(as.character(linking.key$QTL_No))

ld.SNPs = merge(ld.SNPs, linking.key, by = c("PeakChr", "PeakPos", "PeakAllele"))
ld.SNPs$PeakLoc = factor(as.character(paste(ld.SNPs$PeakChr, ld.SNPs$PeakPos, sep = ".")))

GWAS.hits = merge(GWAS.hits, linking.key, by = c("QTL_No"), all.x = TRUE)
GWAS.hits = merge(GWAS.hits, subset(ld.SNPs, , select = c("Chr", "Pos", "QTL", "PeakChr", "PeakPos", "Rsq")), by.x = c("Chr", "bp", "QTL", "PeakChr", "PeakPos"), by.y = c("Chr", "Pos", "QTL", "PeakChr", "PeakPos"), all.x = TRUE)
GWAS.hits$PeakLoc = factor(as.character(paste(GWAS.hits$PeakChr,GWAS.hits$PeakPos, sep = ".")))
GWAS.hits$Trait = factor(GWAS.hits$Trait, levels = trait.order)

genes.coord = merge(genes.coord, linking.key, by = c("GRMZM_ID", "QTL"))

expr.corr$QTL_No = factor(as.character(expr.corr$QTL_No), levels = levels(linking.key$QTL_No))
expr.corr$root = NULL
expr.corr$shoot = NULL
expr.corr.melt = melt(expr.corr, id.vars = c("QTL_No","QTL", "GRMZM_ID", "Trait"), variable.name = "Sampling", value.name = "Corr")
expr.corr.pval$QTL_No = factor(as.character(expr.corr.pval$QTL_No), , levels = levels(linking.key$QTL_No))
expr.corr.pval.melt = melt(expr.corr.pval, id.vars = c("QTL_No", "QTL", "GRMZM_ID", "Trait"), variable.name = "Sampling", value.name = "CorrPval")
QTL.expr = merge(expr.corr.melt, expr.corr.pval.melt, by = c("QTL_No", "QTL", "GRMZM_ID", "Trait", "Sampling"))
QTL.expr = merge(QTL.expr, linking.key, by = c("QTL_No", "QTL", "GRMZM_ID")) 
QTL.expr$Sampling = as.character( sub("^X", "", QTL.expr$Sampling))
QTL.expr$Trait = factor(QTL.expr$Trait, levels = trait.order)
QTL.expr$Sampling = factor(QTL.expr$Sampling, levels = samples.order)
QTL.expr$SamplingLabels = sub("DAP", "", QTL.expr$Sampling)


# Set GWAS and QTL expression trait formats
GWAS.hits = merge(GWAS.hits, trait.format, by = c("Trait"))
QTL.expr = merge(QTL.expr, trait.format, by = c("Trait"))


# Set genes direction dash positions centered on the middle of the gene
genes.coord$DashXStart = (genes.coord$start + genes.coord$stop)/2 + as.integer(paste(genes.coord$strand, ceiling(kGeneArrowProp * (genes.coord$stop - genes.coord$start)), sep = ""))
genes.coord$DashXEnd = (genes.coord$start + genes.coord$stop)/2 - as.integer(paste(genes.coord$strand, ceiling(kGeneArrowProp * (genes.coord$stop - genes.coord$start)), sep = ""))

# Select SNP subsets within the given window from the peak
ldSNPs.window = subset(ld.SNPs, Chr == PeakChr & Pos >= (PeakPos - kSearchSpace) & Pos <= (PeakPos + kSearchSpace) & Rsq >= 0)
GWAS.window = subset(GWAS.hits, Chr == PeakChr & bp >= (PeakPos - kSearchSpace) & bp <= (PeakPos + kSearchSpace))



# Bin the ld window for ld smoothing and SNP density calc
ldSNPs.window.bins = data.frame()

for (peak.pos in as.vector(unique(ldSNPs.window$PeakLoc))){
    peak.ldSNPs = subset(ldSNPs.window, PeakLoc == peak.pos);
    nBreaks = 1 + ceiling((max(peak.ldSNPs$Pos) - min(peak.ldSNPs$Pos)) / kBinSize)
    peak.ldSNPs$Breaks = cut(peak.ldSNPs$Pos, nBreaks, include.lowest = FALSE);

    nr.snps = nrow(peak.ldSNPs);
    
    peak.ldSNPs.bins = unique(subset(peak.ldSNPs, , select = c("Breaks", "PeakChr", "PeakPos", "QTL_No", "QTL", "GRMZM_ID", "PeakLoc")))
    peak.ldSNPs.bins$lower = as.numeric(sub("\\((.+),.*", "\\1", peak.ldSNPs.bins$Breaks))
    peak.ldSNPs.bins$upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", peak.ldSNPs.bins$Breaks))
    peak.ldSNPs.bins$midPos = (peak.ldSNPs.bins$lower + peak.ldSNPs.bins$upper)/2
    
    peak.ldSNPs.bins = merge(peak.ldSNPs.bins, setNames(aggregate(data = peak.ldSNPs, Rsq ~ Breaks, mean, na.rm = FALSE), c("Breaks", "meanRsq")), by = c("Breaks"), all.x = TRUE);
    peak.ldSNPs.bins = merge(peak.ldSNPs.bins, setNames(aggregate(data = peak.ldSNPs, Pos ~ Breaks, FUN = function(x){sum( !is.na(x) )}), c("Breaks", "SNPcount")), by = c("Breaks"),  all.x = TRUE);
    peak.ldSNPs.bins = merge(peak.ldSNPs.bins, setNames(aggregate(data = peak.ldSNPs, Pos ~ Breaks, min, na.rm = TRUE), c("Breaks", "minPos")), by = c("Breaks"),  all.x = TRUE);
    peak.ldSNPs.bins = merge(peak.ldSNPs.bins, setNames(aggregate(data = peak.ldSNPs, Pos ~ Breaks, max, na.rm = TRUE), c("Breaks", "maxPos")), by = c("Breaks"),  all.x = TRUE);
    peak.ldSNPs.bins = merge(peak.ldSNPs.bins, setNames(aggregate(data = peak.ldSNPs, Rsq ~ Breaks, max, na.rm = TRUE), c("Breaks", "maxRsq")), by = c("Breaks"),  all.x = TRUE);

    ldSNPs.window.bins = rbind(ldSNPs.window.bins, peak.ldSNPs.bins);
}


# find contiguous regions of no SNPs
#binsAllBreaks=setNames(as.data.frame(levels(ldSNPs.window.bins$Breaks)), c("Breaks"))
#binsAllBreaks$lower=as.numeric( sub("\\((.+),.*", "\\1", binsAllBreaks$Breaks) )
#binsAllBreaks$upper=as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", binsAllBreaks$Breaks) )
#binsMissing=subset(binsAllBreaks, !is.element(Breaks, as.vector(unique(ldSNPs.window.bins$Breaks))))
#binsMissing.melt=melt(binsMissing, id.vars=c("Breaks"), measure.vars=c("lower","upper"));
#binsGapEdges=subset(merge(binsMissing.melt, subset(aggregate(variable ~ value, data = binsMissing.melt, length), variable == 1), by=c("value")),,select=c("value","variable.x"))
#binsGaps=cbind(setNames(subset(binsGapEdges, variable.x == "lower", select=c("value")), c("lower")), setNames(subset(binsGapEdges, variable.x == "upper", select=c("value")), c("upper")))





##################################
##################################
#
# PLOTTING FUNCTIONS
#
##################################
##################################



#=================================
#=================================
#
# Gene expression plot
#
#=================================
#=================================

exprPlot = function(expression.toPlot, GWAS.traits){

    # initialize subplots
    exprSignifPlot  = ggplot()
    exprPlot = ggplot()
    exprLabelsPlot = ggplot()

    # plot only if we have relevant data
    if( nrow(expression.toPlot) ){

        sampling.labels = unique(subset(expression.toPlot, , select = c("Sampling", "SamplingLabels")))
        sampling.labels.list = setNames(as.character(sampling.labels$SamplingLabels), sampling.labels$Sampling)
        
        expr.xlim.corr = c(min(as.integer(expression.toPlot$Sampling)), max(as.integer(expression.toPlot$Sampling)))
        expr.xlim.labels = c(0, 0.5)
        
        expr.ylim.corr = c(-1, 1)
        expr.ylim.sig = c(1, length(trait.group.order) + 1)  # default, will be adjusted from data below
        expr.ylim.labels = expr.ylim.corr

        expr.nrows = length(unique(expression.toPlot$Trait))
        expr.point.size = 3

        # Create different labelling aesthetics depending on whether the trait is in the GWAS plot or not
        #-------------------------------------------------------------------------------------------------
        aes.col = data.frame(
            "LabelText" = rep(NA, expr.nrows),
            "LabelColour" = rep(NA, expr.nrows),
            "LabelFill" = rep(NA, expr.nrows),
            "LabelAlpha" = rep(NA, expr.nrows),
            "PointColour" = rep(NA, expr.nrows),
            "PointFill" = rep(NA, expr.nrows),
            "PointShape" = rep(NA, expr.nrows),
            "PointSize" = rep(NA, expr.nrows),
            "PointAlpha" = rep(NA, expr.nrows),
            "PointStrokeSize" = rep(NA, expr.nrows),
            "LineColour" = rep(NA, expr.nrows),
            "LineType" = rep(NA, expr.nrows),
            "LineSize" = rep(NA, expr.nrows),
            "LineAlpha" = rep(NA, expr.nrows)
            )
        
        expr.toPlot = cbind(expression.toPlot, aes.col)


        # Traits also present in the ld plot
        #-----------------------------------
        expr.toPlot.LdPlot = subset(expr.toPlot, is.element(Trait, GWAS.traits))
        if(nrow(expr.toPlot.LdPlot)){
        # label aesthetics
            expr.toPlot.LdPlot$LabelText = paste("plain(", expr.toPlot.LdPlot$TraitLabel, ")", sep="")
            expr.toPlot.LdPlot$LabelColour = "#000000"
            expr.toPlot.LdPlot$LabelFill = "#FFFFFF"
            expr.toPlot.LdPlot$LabelAlpha = as.numeric(1)
        # point aesthetics
            expr.toPlot.LdPlot$PointColour = "#000000"
            expr.toPlot.LdPlot$PointFill = expr.toPlot.LdPlot$TraitColour
            expr.toPlot.LdPlot$PointShape = expr.toPlot.LdPlot$TraitShape
            expr.toPlot.LdPlot$PointSize = as.numeric(2)
            expr.toPlot.LdPlot$PointAlpha = as.numeric(1)
            expr.toPlot.LdPlot$PointStrokeSize = as.numeric(0.3)
        # line aesthetics
            expr.toPlot.LdPlot$LineColour = expr.toPlot.LdPlot$TraitColour
            expr.toPlot.LdPlot$LineType = "solid"
            expr.toPlot.LdPlot$LineSize = as.numeric(0.5)
            expr.toPlot.LdPlot$LineAlpha = as.numeric(1)
        }

       

        # Traits not present in the LD plot
        #----------------------------------
        expr.toPlot.noLdPlot = subset(expression.toPlot, ! is.element(Trait, GWAS.traits))
        if(nrow(expr.toPlot.noLdPlot)){
        # label aesthetics
            expr.toPlot.noLdPlot$LabelText = paste("italic(", expr.toPlot.noLdPlot$TraitLabel, ")", sep="")
            expr.toPlot.noLdPlot$LabelColour = "#888888"
            expr.toPlot.noLdPlot$LabelFill = "#FFFFFF"
            expr.toPlot.noLdPlot$LabelAlpha = as.numeric(1)
        # points aesthetics
            expr.toPlot.noLdPlot$PointColour = expr.toPlot.noLdPlot$TraitColour
            expr.toPlot.noLdPlot$PointFill = "#FFFFFF"
            expr.toPlot.noLdPlot$PointShape = expr.toPlot.noLdPlot$TraitShape
            expr.toPlot.noLdPlot$PointSize = as.numeric(2)
            expr.toPlot.noLdPlot$PointAlpha = as.numeric(1)
            expr.toPlot.noLdPlot$PointStrokeSize = as.numeric(0.5)
        # line aesthetics
            expr.toPlot.noLdPlot$LineColour = expr.toPlot.noLdPlot$TraitColour
            expr.toPlot.noLdPlot$LineType = "dashed"
            expr.toPlot.noLdPlot$LineSize = as.numeric(0.5)
            expr.toPlot.noLdPlot$LineAlpha = as.numeric(1)
        }

        # combine the subsets
        expr.toPlot = rbind(expr.toPlot.noLdPlot, expr.toPlot.LdPlot)


        ########################
        # Create the main plots
        ########################


        #===============================
        # Signif. expr. labels sub-plot
        #===============================

        expr.corr.signif = subset(expr.toPlot, CorrPval <= kExprPvalCutoff)

        # If we have significant correlation, group traits by trait group and create labels
        #----------------------------------------------------------------------------------
        if (nrow(expr.corr.signif)){            
            expr.corr.signif.labels = aggregate(
                TraitLabel ~ Sampling + QTL + TraitGroup,
                data = expr.corr.signif,
                FUN = function(x){
                    paste(
                        "list(",
                        paste(x, collapse = ", "),
                        ")",
                        sep="")
                }
                )
            
            expr.corr.signif.labels$Sampling = factor(expr.corr.signif.labels$Sampling, levels=samples.order)
            expr.corr.signif.labels$TraitGroup = factor(expr.corr.signif.labels$TraitGroup, levels=trait.group.order)


            # change the proportion of the significance plot based on the nr. of rows to be plotted
            #--------------------------------------------------------------------------------------
            max.signif.nrows = length(trait.group.order)
            signif.nrows = length(as.vector(unique(expr.corr.signif.labels$TraitGroup)))
            max.signif.prop = expr.plot.prop.row[1]
            signif.prop = (max.signif.prop / max.signif.nrows) * signif.nrows
            expr.plot.prop.row = c(signif.prop, 1 - signif.prop)


            # plot the labels
            exprSignifPlot = exprSignifPlot +
                geom_text(
                    data = expr.corr.signif.labels,
                    aes(
                        x = as.integer(Sampling),
                        y = TraitGroup,
                        label = TraitLabel
                        ),
                    parse = TRUE,
                    size = 2,
                    hjust = 0.5,
                    vjust = 0
                    )
        }
        else{ #if no significant hits make the correlation plot full height
            expr.plot.prop.row = c(0, 1)
            signif.nrows = 0
        }

        # Adjust the y axis limits
        expr.ylim.sig = c(1, signif.nrows + 1)
        
        
    #================
    # Expr. sub-plot
    #================

        exprPlot = exprPlot +

            # horizontal line for minimum abs(r): r = 0
            geom_hline(yintercept = 0, linetype = "longdash", size = 0.2, colour = "gray80") +

            # horizontal lines for maximum abs(r): r = +/- 1
            geom_hline(yintercept = c(-1, 1), linetype = "solid", size = 0.1, colour = "gray80") +

            # trend lines linking the ordered sampling time points
            geom_line(
                data = expr.toPlot,
                aes(
                    group = Trait,
                    x = Sampling,
                    y = Corr,
                    colour = as.character(LineColour),
                    size = as.numeric(LineSize),
                    linetype = as.character(LineType),
                    alpha = as.numeric(LineAlpha)
                    )
                ) +

            # symbol points for specific Trait X Sampling r values
            geom_point(
                data = expr.toPlot,
                aes(
                    x = Sampling,
                    y = Corr,
                    shape = as.integer(PointShape),
                    colour = as.character(PointColour),
                    size = as.numeric(PointSize),
                    stroke = as.numeric(PointStrokeSize),
                    fill = as.character(PointFill)                    
                    )
                )

       

    
        #=========================
        # Expression trait labels
        #=========================

        # set the height of the label anchor based on the rightmost sampling point
        #-------------------------------------------------------------------------
        expr.label.toPlot = subset(expr.toPlot, as.integer(Sampling) == expr.xlim.corr[2])

        exprLabelsPlot = exprLabelsPlot +
            
            geom_text_repel(
                data = expr.label.toPlot,
                    x = 0,
                    aes(
                        y = Corr,
                        label = as.character(LabelText),
                        fill = as.character(LabelFill),
                        colour = as.character(LabelColour),
                        alpha = as.numeric(LabelAlpha)
                        ),
                    nudge_x = 0.25,
                    segment.size = 0.25,
                    size = 2,
                    box.padding = unit(0.1, "lines"),
                    point.padding = unit(0, "lines"),
                    segment.color="black",
                    parse = TRUE
                    )
    }
    


    ##########################################
    # Set overall formatting for the subplots
    ##########################################

    
    #========================================
    # Format expression significance subplot
    #========================================
    
    exprSignifPlot = exprSignifPlot +
        coord_cartesian(xlim = expr.xlim.corr, ylim = expr.ylim.sig) +
        theme_tufte(ticks = FALSE) +
        scale_linetype_identity() +
        scale_fill_identity() +
        scale_colour_identity() +
        scale_shape_identity() +
        scale_size_identity() +
        scale_x_discrete(expand = c(0, 0.1)) +  # expand sligthly to fit x axis labels
        scale_y_discrete(expand = c(0, 0)) +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
#            axis.title.y = element_text(size = kAxisLabelSize, colour = NA, face = kAxisLabelFace),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0,0,0,0), "line")
            ) +
        labs(title = NULL, x = NULL, y = "TraitGroup")


    #===========================
    # Format expression subplot
    #===========================
    
    exprPlot = exprPlot +
        coord_cartesian(xlim = expr.xlim.corr, ylim = expr.ylim.corr) +
        theme_tufte(ticks = FALSE) +
        scale_linetype_identity() +
        scale_fill_identity() +
        scale_colour_identity() +
        scale_shape_identity() +
        scale_size_identity() +
        scale_x_discrete(expand = c(0,0.1), labels = sampling.labels.list) +  # expand slightly to fit x axis labels
        scale_y_continuous(expand = c(0,0)) +
        theme(
            legend.position = "none",
#                axis.text.x = element_blank(),
            axis.text.x = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
            axis.title.x = element_blank(),
#                axis.title.x = element_text(size = kAxisLabelSize, colour = kAxisLabelColour, face = kAxisLabelFace),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0,0,0,0), "line")
            ) +
        labs(title = NULL, x = NULL, y="r(FPKM, JL AEE)")

    
    #==================================
    # Format expression labels subplot
    #==================================
    
    exprLabelsPlot = exprLabelsPlot +
        coord_cartesian(xlim = expr.xlim.labels, ylim = expr.ylim.labels) +
        theme_tufte(ticks = FALSE) +
        scale_linetype_identity() +
        scale_fill_identity() +
        scale_colour_identity() +
        scale_shape_identity() +
        scale_size_identity() +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme(
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = kAxisTextSize, colour = NA, face = kAxisTextFace),
            axis.title.x = element_text(size = kAxisLabelSize, colour = NA, face = kAxisLabelFace),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0,0,0,0), "line")
            ) +
        labs(title = NULL, x = "Traits", y = NULL)
        

    
    ########################
    # Set x axis label text
    ########################

    xlab.grob = textGrob(
        label = "Days after pollination (DAP)",
        hjust = 0.5,
        vjust = 0,
        just = "centre",
        gp = gpar(fontsize = kAxisLabelSize, colour = kAxisLabelColour, face = kAxisLabelFace)
        )


    
    ########################
    # Combine the sub-plots
    ########################
    
    plot = arrangeGrob(
        
        # sub-plots, in grid layout
        exprSignifPlot, plot_blank(),
        exprPlot, exprLabelsPlot,

        # formatting parameters
        padding = unit(1,"line"),
        ncol = 2,
        nrow = 2,
        widths = expr.plot.prop.col,
        heights = expr.plot.prop.row,
        bottom = xlab.grob
        )

    ##################################################
    # Uncomment below to add a border around the plot
    # plot = gtable::gtable_add_grob(plot, rectGrob(gp=gpar(lwd=kBorderSize, fill=NA, colour=kBorderColour)), 1, 2, 3, 1)
    ##################################################

    # return a grob element
    return(plot)
}







#=================================
#=================================
#
# Main LD plot
#
#=================================
#=================================


ldPlot = function(ldSNPs.toPlot, ldSNPs.toPlot.bins, GWAS.toPlot, genes.toPlot){
    
    
    # Initialize the ld plot
    RMIPPlot = ggplot()
    GenesPlot = ggplot()
    ldGWASPlot = ggplot()

    # Set up plot limits
    plot.xlim = c(min(ldSNPs.toPlot.bins$lower), max(ldSNPs.toPlot.bins$upper))
    plot.ylim.RMIP = c(0, 1)
    plot.ylim.genes = c(-1, 1)
    plot.ylim.LD = c(0, 1)

    # Create subsets for the peak markers and r^2 values
    GWAS.Peak.toPlot = aggregate( RMIP ~ PeakChr + PeakPos, data = subset(GWAS.toPlot, Chr == PeakChr & bp == PeakPos ), FUN = max)
    GWAS.Rsq.toPlot = subset(GWAS.toPlot, ! is.na(Rsq))



    
    ########################
    # Create the main plots
    ########################

    
    #===================
    # GWAS RMIP subplot
    #===================
    
    if( nrow(GWAS.toPlot) > 0 ){

        # Set up plot vertical edges
        #----------------------------
        RMIPPlot = RMIPPlot +

            # line for RMIP = 0
            geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "black") +

            # line for RMIP = 1
            geom_hline(yintercept = 1, linetype = "solid", size = 0.1, colour = "gray80")

        
        # Plot the RMIP data
        #--------------------
        RMIPPlot = RMIPPlot +

            # vertical bars with RMIP values
            geom_segment(
                data = GWAS.toPlot,
                y = 0,
                aes(
                    x = bp,
                    xend = bp,
                    yend = RMIP
                    ),
                colour = "black",
                size = 0.3,
                linetype = 1
                )
    

        # Add a custom vertical bar for the peak value
        #-------------------------------------------
        if( nrow(GWAS.Peak.toPlot) > 0 ){
            
            RMIPPlot = RMIPPlot +

                geom_segment(
                    data = GWAS.Peak.toPlot,
                    aes(
                        x = PeakPos,
                        xend = PeakPos,
                        yend = RMIP
                        ),
                    y = 0,
                    colour = "red",
                    size = 0.5
                    ) 
            
        }


        # Add white circles at the base of the vertical bars
        #---------------------------------------------------
        RMIPPlot = RMIPPlot +

            geom_point(
                data = GWAS.toPlot,
                position = "identity",
                aes(
                    x = bp
                    ),
                y = 0,
                size = 1,
                colour="black",
                fill = "white",
                shape = 21,
                stroke = 0.25
                ) 


        # Draw symbols and label the RMIP values
        #----------------------------------------
        RMIPPlot = switch(label.type,

            # non-overlapping, auto-adjusted labels with arrows
            'arrow' = {
                RMIPPlot +
                geom_label_repel(
                    data = GWAS.toPlot,
                    stat = "identity",
                    force = 1,
                    nudge_y = 0.2,
                    nudge_x = 0.2,
                    aes(
                        x = bp,
                        y = RMIP,
                        label = TraitLabel
                        ),
                    colour = "black",
                    fill = "white",
                    alpha = 0.5,
                    segment.color = "grey30",
                    segment.size = 0.2,
                    arrow = arrow(
                        length = unit(0.025, "npc"),
                        ends = "last",
                        type = "closed",
                        ),
                    box.padding = unit(0.3, "lines"),
                    label.padding = unit(0.1, "lines"),
                    point.padding = unit(0.8, "lines"),
                    label.r = unit(0.1, "lines"),
                    label.size = 0,
                    na.rm = TRUE,
                    size = 2.5,
                    fontface = "italic",
                    parse = TRUE
                    ) +
                geom_point(
                    data = GWAS.toPlot,
                    position="identity",
                    aes(
                        x = bp,
                        y = RMIP,
                        fill = as.character(TraitColour),
                        shape = as.integer(TraitShape)
                        ),
                    size = 2,
                    colour="black",
                    stroke=0.3
                    ) 
            },

            # text labels right above the symbol; can/will overlap of plot out of bounds
            'text' = {
                RMIPPlot +
                geom_point(
                    data = GWAS.toPlot,
                    position="identity",
                    aes(
                        x = bp,
                        y = RMIP,
                        fill = as.character(TraitColour),
                        shape = as.integer(TraitShape)
                        ),
                    size = 2,
                    colour="black",
                    stroke=0.3
                    ) +
                geom_label(
                    data = GWAS.toPlot,
                    stat = "identity",
                    nudge_y = 0.09,
                    aes(
                        x = bp,
                        y = RMIP,
                        label = TraitLabel
                        ),
                    colour = "black",
                    fill = "white",
                    alpha = 0.5,
                    label.padding = unit(0.05, "lines"),
                    label.r = unit(0.1, "lines"),
                    label.size = 0,
                    na.rm = TRUE,
                    size = 2.5,
                    fontface = "italic",
                    parse = TRUE
                    )
            },

            # no labels, just symbols
            'off' = {
                RMIPPlot +
                geom_point(
                    data = GWAS.toPlot,
                    position="identity",
                    aes(
                        x = bp,
                        y = RMIP,
                        fill = as.character(TraitColour),
                        shape = as.integer(TraitShape)
                        ),
                    size = 2,
                    colour="black",
                    stroke=0.3
                    )
            }
            )
                   

    }

    
    
    #=====================
    # Gene models subplot
    #=====================

    if(nrow(genes.toPlot)){

        GenesPlot = GenesPlot +

            # the main gene model
            geom_rect(
                data = genes.toPlot,
                aes(
                    xmin = start,
                    xmax = stop,
                    ymin = start * 0 - 1, # need to set this inside aes or x axis won't line up properly
                    ymax = stop * 0 + 1   # need to set this inside aes or x axis won't line up properly
                    ),
                fill = "black",
                colour = NA,
                size = 0.1,
                ) +

            # arrow showing gene model direction
            geom_segment(
                data = genes.toPlot,
                aes(
                    x = DashXEnd,
                    xend = DashXStart,
                    y = start * 0 - 1,   # need to set this inside aes or x axis won't line up properly
                    yend = stop * 0      # need to set this inside aes or x axis won't line up properly
                    ),
                colour = "white",
                size = 0.3,
                ) +
            geom_segment(
                data = genes.toPlot,
                aes(
                    x = DashXStart,
                    xend = DashXEnd,
                    y = start * 0,        # need to set this inside aes or x axis won't line up properly
                    yend = stop * 0 + 1   # need to set this inside aes or x axis won't line up properly
                    ),
                colour = "white",
                size = 0.3,
                ) +

            # border around the gene model to trim the edges of the arrow
            geom_rect(
                data = genes.toPlot,
                aes(
                    xmin = start,
                    xmax = stop,
                    ymin = start * 0 - 1, # need to set this inside aes or x axis won't line up properly
                    ymax = stop * 0 + 1   # need to set this inside aes or x axis won't line up properly
                    ),
                fill = NA,
                colour = "black",
                size = 0.2,
                ) 
    }



    #=======================
    # GWAS LD (r^2) subplot
    #=======================

    #########################################################################
    # Uncomment below to show all LD values (useful for smoothing evaluation)
    # ldGWASPlot = ldGWASPlot + geom_segment(data = ldSNPs.toPlot, y = 0, aes(x = Pos, xend = Pos, yend = Rsq), colour = "grey80", size = 0.1)
    #########################################################################

    
    if( nrow(GWAS.Rsq.toPlot) > 0 ){

        # Set up plot vertical edges
        #----------------------------
        ldGWASPlot = ldGWASPlot +

            # line for r^2 = 0
            geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "black") +

            # lline for r^2 = 1
            geom_hline(yintercept = 1, linetype = "solid", size = 0.1, colour = "gray80")


        # Plot the r^2 data
        #-------------------
        ldGWASPlot = ldGWASPlot +
            
            # Vertical bars with r^2 values
            geom_segment(
                data = GWAS.Rsq.toPlot,
                aes(
                    x = bp,
                    xend = bp,
                    y = Rsq * 0, # need to set this inside aes or x axis won't line up properly
                    yend = Rsq
                    ),
                colour = "black",
                size = 0.3,
                alpha = 1
                ) 
    
        # Add a custom vertical bar for the peak value
        #-------------------------------------------    
        if( nrow(GWAS.Peak.toPlot) > 0 ){
            
            ldGWASPlot = ldGWASPlot +
            
                geom_segment(
                    data = GWAS.Peak.toPlot,
                    aes(
                        x = PeakPos,
                        xend = PeakPos,
                        y = PeakPos * 0,        # need to set this inside aes or x axis won't line up properly
                        yend = PeakPos * 0 + 1  # need to set this inside aes or x axis won't line up properly
                        ),
                    colour = "red",
                    size = 0.5,
                    linetype = "solid"
                    )
        }


        # Add white circles at the base of the vertical bars
        #---------------------------------------------------
        ldGWASPlot = ldGWASPlot +

            geom_point(
                data = GWAS.toPlot,
                position="identity",
                aes(
                    x = bp
                    ),
                y = 0,
                size = 1,
                colour="black",
                fill = "white",
                shape = 21,
                stroke=0.25
                )
    }





    ##########################################
    # Set overall formatting for the subplots
    ##########################################

    
    #========================================
    # Format RMIP sub-plot
    #========================================

    RMIPPlot = RMIPPlot +
        coord_cartesian(ylim = plot.ylim.RMIP, xlim = plot.xlim) +
        theme_tufte(ticks = FALSE) +
        scale_colour_identity() +
        scale_fill_identity() +
        scale_shape_identity() +
        scale_x_continuous(expand=c(0, 0)) +
        scale_y_continuous(expand=c(0, 0.025)) +  # expand slightly so the white circles at the base are not cut off
        theme(
            legend.position="none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
#            axis.text.y = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
            axis.title.y = element_blank(),
#            axis.title.y = element_text(size = kAxisLabelSize, colour = kAxisLabelColour, face = kAxisLabelFace),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0, 0, 0, 0), "line")
            ) +
        labs(title = NULL, x = NULL, y = "RMIP")



    #========================================
    # Format gene models sub-plot
    #========================================

    GenesPlot = GenesPlot +
        coord_cartesian(ylim = plot.ylim.genes, xlim = plot.xlim) +
            theme_tufte(ticks = FALSE) +
            scale_fill_identity() +
            scale_colour_identity() +
            scale_shape_identity() +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                legend.position = "none",
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.y = element_blank(),
#                axis.text.y = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace)
                axis.title.y = element_blank(),
#                axis.title.y = element_text(size = kAxisLabelSize, colour = NA, face = kAxisLabelFace),
                panel.margin = unit(0, "mm"),
                plot.background = element_rect(fill = NA, colour = NA),
                plot.margin = unit(c(0, 0, 0, 0), "line")
            ) +
            labs(title = NULL, x = NULL, y = "Gene")



    #========================================
    # Format GWAS LD sub-plot
    #========================================

    ldGWASPlot = ldGWASPlot +
        coord_cartesian(ylim = plot.ylim.LD, xlim = plot.xlim) +
            theme_tufte(ticks = FALSE) +
            scale_colour_identity() +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_reverse(expand = c(0, 0.025)) + # expand slightly so the white circles at the base are not cut off
            theme(
                legend.position = "none",
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.y = element_blank(),
#                axis.text.y = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
                axis.title.y = element_blank(),
#                axis.title.y = element_text(size = kAxisLabelSize, colour = kAxisLabelColour, face = kAxisLabelFace),
                panel.margin = unit(0, "mm"),
                plot.background = element_rect(fill = NA, colour = NA),
                plot.margin = unit(c(0, 0, 0, 0), "line")
            ) +
            labs(title = NULL, x = NULL, y = expression(r^2*(list(SNP, PeakSNP))))


    
    ########################
    # Combine the sub-plots
    ########################

    plot = arrangeGrob(

        # sub-plots, in grid layout
        RMIPPlot, plot_blank(),
        GenesPlot, plot_blank(),
        ldGWASPlot, plot_blank(),

        # formatting parameters
        ncol = 2,
        nrow = 3,
        widths = ld.plot.prop.col,
        heights = ld.plot.prop.row,
        padding=unit(1,"line")
        )


    ##################################################
    # Uncomment below to add a border around the plot    
    # plot=gtable::gtable_add_grob(plot, rectGrob(gp=gpar(lwd=kBorderSize, fill=NA, colour=kBorderColour)), 1, 2, 3, 1)       
    ##################################################

    # return a grob element
    return(plot);
}






#=================================
#=================================
#
# LD bins and SNP density ribbons
#
#=================================
#=================================


ldRibbonPlot = function(ldSNPs.toPlot.bins){

    # Initialize plot objects
    ldDecayPlot = ggplot()
    coordPlot = ggplot()
    ldSNPsPlot = ggplot()

    # Set up graph aestetic adjustments
    coord.labels.no = 5

    # Calculate limits
    max.SNP = max(ldSNPs.toPlot.bins$SNPcount)
    if (max.SNP < 0) { # if no SNPs
        max.SNP = 0
    }

    plot.xlim = c(min(ldSNPs.toPlot.bins$lower), max(ldSNPs.toPlot.bins$upper))

    plot.ylim.LD = switch(ribbon.type,
        'flat' = c(0, 1),
        'hist' = c(0, 1)
        )

    plot.ylim.SNP = switch(ribbon.type,
        'flat' = c(0, 1),
        'hist' = c(0, max.SNP)
        )

    if(nrow(ldSNPs.toPlot.bins) <= 0){
        plot.prop.row = c(0,1,0)
    }

    plot.ylim.coord = c(0, 0)
    
    # Calculate coordinate labels and breaks, dropping the edge values for aesthetic reasons
    coord = seq(min(ldSNPs.toPlot.bins$lower), max(ldSNPs.toPlot.bins$upper), length.out = coord.labels.no + 2)
    coord = coord[2:(length(coord) - 1)]
    coord.seq = data.frame(
        xPos = coord,
        yPos = rep(0, length(coord))
        )
    
   

    
    if( nrow(ldSNPs.toPlot.bins) ){

    #====================
    # Binned LD sub-plot
    #====================
        
        # Set LD plot based on ribbon.type switch variable
        #--------------------------------------------------
        ldDecayPlot = switch(ribbon.type,

            'flat' = {
                ldDecayPlot +
                geom_rect(
                    data = ldSNPs.toPlot.bins,
                    ymin = 0,
                    aes(
                        xmin = lower,
                        xmax = upper,
                        fill = maxRsq,
                        ymax = 0 * maxRsq + 1 # need to specify in aes as workaround for x axis offset "bug"
                        ),
                    colour=NA
                    )
            },
            
            'hist' = {
                ldDecayPlot +
                geom_hline(yintercept = 1, linetype = "solid", size = 0.1, colour = "gray80") + 
                geom_rect(
                    data = ldSNPs.toPlot.bins,
                    ymin = 0,
                    aes(
                        xmin = lower,
                        xmax = upper,
                        fill = maxRsq,
                        ymax = maxRsq
                        ),
                    colour=NA
                    )
            }
            )
        
                
        
    #=============================
    # Binned SNP density sub-plot
    #=============================
        
        # Set SNP density plot based on ribbon.type switch variable
        #-----------------------------------------------------------
        ldSNPsPlot = switch(ribbon.type,

            'flat' = {
                ldSNPsPlot +
                geom_rect(
                    data = ldSNPs.toPlot.bins,
                    ymin = 0,
                    aes(
                        xmin = lower,
                        xmax = upper,
                        fill = SNPcount,
                        ymax = SNPcount * 0 + 1 # need to specify in aes as workaround for x axis offset "bug"
                        ),
                    colour=NA
                    )
            },
            
            'hist' = {
                ldSNPsPlot +
                geom_hline(
                    yintercept = max.SNP,
                    linetype = "solid",
                    size = 0.1,
                    colour = "gray80"
                    ) +            
                geom_rect(
                    data = ldSNPs.toPlot.bins,
                    ymin = 0,
                    aes(
                        xmin = lower,
                        xmax = upper,
                        fill = SNPcount,
                        ymax = SNPcount
                        ),
                    colour=NA
                    )
            }
            )
    }

    
    #=============================
    # Genome coordinates plot
    #=============================

    coordPlot = coordPlot +
        
        geom_text(
            data = coord.seq,
            aes(
                x = xPos,
                y = yPos,
                label = rescaleMb(xPos)
                ),
            check_overlap = TRUE,
            size = 2,
            colour = kAxisTextColour,
            fontface = kAxisTextFace
            )

    

    ##########################################
    # Set overall formatting for the subplots
    ##########################################

    
    #=========================
    # Format LD bins sub-plot
    #=========================

    ldDecayPlot = ldDecayPlot +
        coord_cartesian(xlim = plot.xlim, ylim = plot.ylim.LD) +
        theme_tufte(ticks = FALSE) +
        scale_fill_distiller(palette = "Greys", limits = c(0, 1), direction = 1, na.value = "white", guide = FALSE) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
#                axis.text.y = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
            axis.title.y = element_blank(),
#                axis.title.y = element_text(size = kAxisLabelSize, colour = kAxisTextColour, face = kAxisLabelFace),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0, 0, 0, 0), "line")
            ) +
        labs(title = NULL, x = NULL, y = "r^2")

    

    #=============================
    # Format SNP density sub-plot
    #=============================
    
    ldSNPsPlot = ldSNPsPlot + 
        coord_cartesian(xlim = plot.xlim, ylim = plot.ylim.SNP) +
        theme_tufte(ticks = FALSE) +
        scale_fill_distiller(palette = "Blues", limits = c(0, max.SNP), direction = 1, na.value = "white", guide = FALSE) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme(
            legend.position="none",
            axis.text.x = element_blank(),
#                axis.text.x = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
#                axis.text.y = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
            axis.title.y = element_blank(),
#                axis.title.y = element_text(size = kAxisLabelSize, colour = kAxisTextColour, face = kAxisLabelFace),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0, 0, 0, 0), "line")
            ) +
        labs(title = NULL, x = NULL, y = paste("# SNPs (0, ", max.SNP, ")", sep=""))



    #====================================
    # Format genome coordinates sub-plot
    #====================================

    coordPlot = coordPlot +
        coord_cartesian(xlim = plot.xlim, ylim = plot.ylim.coord) +
        theme_tufte(ticks = FALSE) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
#                axis.text.y = element_text(size = kAxisTextSize, colour = kAxisTextColour, face = kAxisTextFace),
            axis.title.y = element_blank(),
#                axis.title.y = element_text(size = kAxisLabelSize, colour = NA, face = kAxisLabelFace),
            panel.margin = unit(0, "mm"),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = unit(c(0, 0, 0, 0), "line")
            ) +
        labs(title = NULL, x = NULL, y = NULL)


    
    ########################
    # Set x axis label text
    ########################
    
    xlab.grob = textGrob(
        label = paste("Position on chromosome ", as.vector(unique(ldSNPs.toPlot.bins$PeakChr))[1], " (Mb)", sep=""),
        hjust = 0.5,
        vjust = 0,
        just = "centre",
        gp = gpar(fontsize = kAxisLabelSize, colour = kAxisLabelColour, face = kAxisLabelFace)
        )


    
    ########################
    # Combine the sub-plots
    ########################    

    plot = arrangeGrob(

        # sub-plots in grid layout
        ldDecayPlot, plot_blank(),
        coordPlot, plot_blank(),
        ldSNPsPlot, plot_blank(),

        # formatting parameters
        ncol = 2,
        nrow = 3,
        widths = ribbon.plot.prop.col,
        heights = ribbon.plot.prop.row,
        padding=unit(1,"line"),
        bottom = xlab.grob
        )


    ##################################################
    # Uncomment below to add a border around the plot    
    # plot=gtable::gtable_add_grob(plot, rectGrob(gp=gpar(lwd=kBorderSize, fill=NA, colour=kBorderColour)), 1, 2, 4, 1)
    ##################################################

    #return a grob object
    return(plot);

}
    






##################################
##################################
#
# PLOT THE DATA
#
##################################
##################################



# Open output file
#------------------
Cairo(file = out.file, height = 5.6, width = 9.1, units = 'in',type = out.type)


# Loop through the QTLs
#-----------------------
for (QTL.select in as.vector(levels(linking.key$QTL)) ){

    # Select data subset
    #--------------------
    QTL.expr.subset = subset(QTL.expr, QTL == QTL.select)
    ldSNPs.window.subset = subset(ldSNPs.window, QTL == QTL.select)
    ldSNPs.window.bins.subset = subset(ldSNPs.window.bins, QTL == QTL.select)
    GWAS.window.subset = subset(GWAS.window, QTL == QTL.select)
    GWAS.traits = as.vector(unique(GWAS.window.subset$Trait))
    genes.coord.subset = subset(genes.coord, QTL == QTL.select)

    # Set up subplot ratios
    #-----------------------

    # top level plots
    plot.prop.col = c(1)

    plot.prop.row = switch(ribbon.type,
        'flat' = c(0.035, 0.30, 0.02, 0.035, 0.50, 0.11),
        'hist' = c(0.035, 0.30, 0.02, 0.035, 0.40, 0.21)
        )

    if(nrow(ldSNPs.window.bins.subset) <= 0){ # no SNPs ribbon
        plot.prop.row = c(0.035, 0.30, 0.02, 0.035, 0.55, 0.6)
    }


    # expression subplot defaults
    expr.plot.prop.col = c(0.92, 0.08)
    expr.plot.prop.row = c(0.35, 0.65) # will be rescaled based on nr. of signif. trait groups

    # RMIP and LD subplot defaults
    ld.plot.prop.col = c(0.92, 0.08)
    ld.plot.prop.row = c(0.47, 0.06, 0.47)
    
    # ribbons subplot subplot
    ribbon.plot.prop.col = c(0.92, 0.08)    
    ribbon.plot.prop.row = switch(ribbon.type,
        'flat' = c(0.35, 0.30, 0.35),
        'hist' = c(0.43, 0.14, 0.43)
        )



    ####################################
    # Subfig "A" title
    ####################################
    subfigA.title.grob = textGrob(
        label = "A",
        x = unit(0.1, "lines"),
        y = unit(0.1, "lines"),
        hjust = 0,
        vjust = 0,
        gp = gpar(fontsize = kSubTitleSize, colour = kSubTitleColour, face = kSubTitleFace)
        )


    ####################################
    # Subfig "B" title
    ####################################
    subfigB.title.grob = textGrob(
        label = "B",
        x = unit(0.1, "lines"),
        y = unit(0.1, "lines"),
        hjust = 0,
        vjust = 0,
        gp = gpar(fontsize = kSubTitleSize, colour = kSubTitleColour, face = kSubTitleFace)
        )


    ####################
    # Create the figure
    ####################
    
    grid.arrange(

        subfigA.title.grob,
        
        exprPlot(
            QTL.expr.subset,
            GWAS.traits
            ), 
        
        plot_blank(), # separator between A and B sub-figures

        subfigB.title.grob,
        
        ldPlot(
            ldSNPs.window.subset,
            ldSNPs.window.bins.subset,
            GWAS.window.subset,
            genes.coord.subset
            ),

        ldRibbonPlot(
            ldSNPs.window.bins.subset
        ),

        ncol = 1,
        nrow = 6,
        widths = plot.prop.col,
        heights = plot.prop.row,
        top = paste("QTL ", QTL.select, sep="")

        )
        
}

# Close the plotting device (file)
#----------------------------------
dev.off()
