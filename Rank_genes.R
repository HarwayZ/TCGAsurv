library(survival) # Support Kaplan-Meier methods
library(data.table) # Support fread
source('functions.R')

datadir  <- 'Formatted'
cancers <- dir( datadir )
# cancers <- cancers[ -(1:(which(cancers == 'THYM')-1)) ]
cancers <- 'GBM'

for(mycancer in cancers) {

    survfile <- paste(datadir, mycancer, 'survival.txt', sep = '/')
    esetfile <- paste(datadir, mycancer, 'expression.txt', sep = '/')
    clinfile <- paste(datadir, mycancer, 'clinical.txt', sep = '/')
    clin2file<- paste(datadir, mycancer, 'Published version of clinical data.csv', sep='/')
    outfile <- paste0('Ranked/',mycancer,'.conf2.txt')

    eset <- fread(esetfile, sep = "\t", header = TRUE, showProgress = FALSE)
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
                       na.strings = c('[Not Evaluated]', '[Not Available]', '[Unknown]', '[Discrepancy]', '[Not Applicable]', '[Completed]'))

    clin2 <- read.csv2(file = clin2file, header = TRUE,
                       na.strings = c('NA', '#I/T', ''))
    clin2$X <- NULL
    clin2$Case.ID <- as.character(clin2$Case.ID)

    #### Remove patients with no survival data
    clin  <-  clin[ clin$barcode %in% surv$Barcode, ]
    clin2 <- clin2[clin2$Case.ID %in% surv$Barcode, ]

    #### Remove patients with invalid serial times
    clin  <-  clin[surv$Time > 0, ]
    surv  <-  surv[surv$Time > 0, ]

    genes <- names(eset)[-1]
    pvals <- numeric(length(genes))
    adjpv <- numeric(length(genes))
    confs <- character(length(genes))
    confs2<- character(length(genes))
    warns <- character(length(genes))
    bests <- character(length(genes))
    nlow  <- numeric(length(genes))
    nhigh <- numeric(length(genes))

    for(i in 1:length(genes)) {

        gene <- genes[i]
        blank <- "                                           "
        if(i == 1) {
            cat(mycancer,": ",i," of ",length(genes),": ",gene,"\r",sep="")
        } else {
            cat(mycancer,": ",i," of ",length(genes),": ",gene,blank,"\r",sep="")
        }
        flush.console()

        #### Separate cohort based on expression of gene ####
        group <- splitdata(eset, gene, surv$Barcode)

        nlow[i]  <- sum(group == "LOW" )
        nhigh[i] <- sum(group == "HIGH")

        if( nhigh[i] < 10 ) {
            warns[i] <- paste("WARNING: Only", nhigh[i],"patients with high gene expression")
        } else if (nlow[i] < 10) {
            warns[i] <- paste("WARNING: Only", nlow[i], "patients with low gene expression")
        } else if( nhigh[i]/nlow[i] > 4 || nlow[i]/nhigh[i] > 4 ) {
            warns[i] <- paste("WARNING: Uneven distribution of patients")
        } else {
            warns[i] <- 'None'
        }

        # Skip plot if all groups are empty
        if( any(c(nlow[i],nhigh[i]) < 2) ) {
            pvals[i] <- NA
            confs[i] <- '-'
           confs2[i] <- '-'
            bests[i] <- '-'
            next
        }

        #### Add column "group" to survival file ####
        surv$Group <- group
        clin$Group <- group
       clin2$Group <- group[surv$Barcode %in% clin2$Case.ID]

        # Remove subjects with no expression group
        subsurv <- subset(surv, Group != 0 )

        #### Calculate difference between groups ?survdiff ####
        temp <- survdiff(Surv(Time, Status) ~ Group, data = subsurv)

        #### Calculate p value from chisq ####
        pval <- 1 - pchisq(temp$chisq, length(temp$n) - 1)
        pval <- signif(pval, digits = 6)
        pvals[i] <- pval

        if( pval <= 0.05 ) {
            ### Test for confounders ###
            conf  <- confounders( clin, gene)
            conf2 <- confounders(clin2, gene)

            conf <- paste(conf, collapse = ";\t")
            conf2 <- paste(conf2, collapse = ";\t")
            if( conf == "" ) {
                confs[i] <- 'None'
            } else {
                confs[i] <- conf
            }
            if( conf2 == "" ) {
                confs2[i] <- 'None'
            } else {
                confs2[i] <- conf2
            }
        }
        if(confs[i] == "") {
            confs[i] <- '-'
        }
        if(confs2[i] == "") {
            confs2[i] <- '-'
        }
        if(temp$obs[1] < temp$exp[1]) {
            bests[i] <- 'High'
        } else if(temp$obs[2] < temp$exp[2]) {
            bests[i] <- 'Low'
        } else {
            print(pval)
            bests[i] <- 'ERROR'
            #### Plot survfit element ?plot ####
            mysurv <- survfit(Surv(Time, Status) ~ Group, data = subsurv)
            mybeginplot(mysurv); myplot(mysurv); addlegend("topright")
            browser()
        }

    }
    cat(mycancer,": ","Writing table...                              \r",sep="")
    flush.console()

    adjpv <- p.adjust(pvals, method = "BH", n = length(na.omit(pvals)))
    adjpv <- signif(adjpv, digits = 6)
    ii <- order(pvals)
    data <- data.frame(genes[ii], bests[ii], pvals[ii], adjpv[ii], nlow[ii], nhigh[ii], warns[ii], confs2[ii])
    colnames(data) <- c('Gene', 'Best surv. group', 'P-value', 'Adj.pval', 'Nlow', 'Nhigh', 'Warnings', 'Extra confounders')

    write.table(data, file = outfile, quote=F, sep="\t", na='-', row.names=F)
}
