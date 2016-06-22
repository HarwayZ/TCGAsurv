##################################
# surv_RNAseq for R              #
# Created by Mathias Husted Torp #
# 24/06-2016                     #
##################################

# Choose cancer types and genes
mycancers    <- FALSE        # Vector of cancers to perform the analysis on. FALSE = ALL
mygenes      <- c('ASS1', 'ARG1', 'ARG2', 'SLC7A1', 'SLC7A2')       # Required. Vector of genes to perform the analysis on. E.g. c('TP53', 'TP53,ASS1')

options(echo=T)
library(survival) # Support Kaplan-Meier methods
library(data.table, quietly = TRUE) # Support fread
source("functions.R")

# Other arguments
trimxaxs     <- FALSE        # Trim the x-axis of each plot, or keep the same length for all
legendpos    <- "bottomleft" # Position of the legend box. ("bottomleft" or "topright")
forceunit    <- "years"        # Force a time unit ("days", "months", "years")
pairwise     <- TRUE        # Do pairwise analysis of all combinations of the input genes
sdweight     <- "none"       # ("none", "half", "full", number) Omit patients with expression close to the median
signif       <- 0.05          # Cases below this P-value are plotted. Set to 1 to plot all. Default: 0.05

# Folder structure
datadir      <- 'Formatted'
outdir       <- 'Output'

if(sdweight == 'None') { sdweight <- 0
} else if ( sdweight == 'Half' ) { sdweight <- 0.5
} else if ( sdweight == 'Full' ) { sdweight <- 1.0
}

if( mygenes == '' || !exists('mygenes') ) {
    stop('No genes specified')
}

# Get list of all genes in the analysis
uniquegenes <- unique(unlist(strsplit(mygenes, ',', fixed = TRUE)))

if( isFALSE(mycancers) ) {
    allcancers <- dir( datadir )
    mycancers = allcancers
    outfile = paste0(outdir, '/All_', paste0(uniquegenes, collapse = '-'), '.pdf' )
} else {
    outfile <- paste0(outdir, '/',
                      paste0(mycancers, '_', collapse = '_'),
                      paste0(uniquegenes, collapse = '-'), '.pdf' )
}

if( isTRUE(pairwise) ) {
    combinations <- combn(x = uniquegenes, m = 2, paste, collapse = ',')
    mygenes <- append(uniquegenes, values = combinations, after = length(uniquegenes))
    ## Append '_pairwise' to outfile name
    outfile <- gsub('^([^[:space:]]+)(\\.pdf)$', '\\1_pairwise\\2', outfile)
}

pdf(file = outfile, width=11, height=8.5, pointsize=12, paper='special')
for(cancer in mycancers) {
    cat(paste("Cancer:", cancer, "\r\n"))

    survfile <- paste(datadir, cancer, 'survival.txt',   sep = '/')
    esetfile <- paste(datadir, cancer, 'expression.txt', sep = '/')
    clinfile <- paste(datadir, cancer, 'clinical.txt',   sep = '/')

    #### Read eset file ####
    eset <- fread(input = esetfile, sep = "\t", header = TRUE, select = c("TCGA_Barcode", uniquegenes), showProgress = FALSE )

    if( length(uniquegenes) != length(eset) - 1 ) {
        stop(paste("The gene", uniquegenes[!uniquegenes %in% names(eset)], "is not in the data set"))
    }

    # Remove outliers based on PCA
    eset <- removeoutliers('Outliers.txt', eset)

    # Extract names and convert to survival-type
    name <- substr(eset$TCGA_Barcode, 0, 12)

    npatients <- length(name)

    #### Read survival file ####
    surv <- read.table(file = survfile, header = TRUE, sep = "\t", row.names = NULL,
                       colClasses = c('character', 'integer', 'integer'),
                       col.names = c('Barcode', 'Time', 'Status')
    )

    #### Read clinical file ####
    clin <- read.table(file = clinfile, header = TRUE, sep = "\t",
                       na.strings = c('[Not Evaluated]', '[Not Available]', '[Unknown]', '[Discrepancy]', '[Not Applicable]'), quote = "")

    #### Remove patients with no survival data
    clin <- clin[clin$barcode %in% surv$Barcode, ]

    #### Remove patients with invalid serial times
    clin <- clin[surv$Time > 0, ]
    surv <- surv[surv$Time > 0, ]

#     if( isTRUE(pairwise) ) {
#         genes <- append(mygenes, after=length(mygenes), combn(mygenes, m = 2, FUN = paste, collapse=','))
#     }

    for(i in 1:length(mygenes)) {
        # Check for paired test
        gene <- unlist(strsplit(mygenes[i], ',', fixed = TRUE))

        #### Separate cohort based on expression of gene ####
        group <- splitdata(eset, gene, surv$Barcode)

        if(length(gene) == 1) {
            nlow  <- sum(group == "LOW" )
            nhigh <- sum(group == "HIGH")

            # Skip plot if all groups are empty
            try(checkforemptygroups(nlow, nhigh, cnc = cancer), next)
        } else if(length(gene) == 2) {
            nlowlow   <- sum(group == "LOWLOW"  )
            nlowhigh  <- sum(group == "LOWHIGH" )
            nhighlow  <- sum(group == "HIGHLOW" )
            nhighhigh <- sum(group == "HIGHHIGH")

            # Print empty plot, if all groups are empty
            try(checkforemptygroups(nlowlow, nlowhigh, nhighlow, nhighhigh), next)
        } else {

        }

        #### Add column "group" to survival file ####
        surv$Group <- group
        clin$Group <- group

        # Remove subjects with no expression group
        subsurv <- subset(surv, Group != 0 )

        #### Change time unit ####
        if(i == 1) {
            tunit <- findtimeunit(tail(subsurv$Time, 1), forceunit)
            if (tunit == 'years') {
                subsurv$Time <- subsurv$Time / 365
                maxtime <- tail(subsurv$Time, 1)
            } else if (tunit == 'months') {
                subsurv$Time <- subsurv$Time / 30.5
                maxtime <- tail(subsurv$Time, 1)
            }
        }

        #### Calculate difference between groups ?survdiff ####
        temp <- survdiff(Surv(Time, Status) ~ Group, data = subsurv)

        #### Calculate p value from chisq ####
        pval <- 1 - pchisq(temp$chisq, length(temp$n) - 1)
        pval <- signif(pval, digits = 6)

        myby <- findmyby(maxtime)
        genelabel <- paste(gene, collapse = ' and ')

        if( pval <= signif ) {
            #### Fit survival ?survfit ####
            mysurv <- survfit(Surv(Time, Status) ~ Group, data = subsurv)

            #### Plot survfit element ?plot ####
            mybeginplot(mysurv, maxtime = maxtime)
            myplot(mysurv, maxtime = maxtime)

            #### Add legend ?legend ?text ####
            addlegend(legendpos)

            ### Test for confounders ###
            # confounders(clin, mygenes[i])
        }
    }
}
invisible(dev.off())
