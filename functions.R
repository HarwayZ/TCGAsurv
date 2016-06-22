findmyby <- function(value) {
    power <- 10^(floor(log10(value)) -1)
    value <- (value / power) * 1.05
    if(value < 20) {
        factor <- 2.5
    } else if (value < 50) {
        factor <- 5
    } else if (value < 100) {
        factor <- 10
    }
    factor * power
}

mybeginplot <- function(mysurv, maxtime = tail(mysurv$time,1), myby = findmyby(maxtime)) {
    mycol <<- c(1, 2, "#21AE21", "#386DFE")
    plot(mysurv, col = 0, axes=F, ylim=c(0,1.05), yaxs='i', xlim=c(0, ceiling(maxtime*1.05)) )
    abline(h=(seq(0.25,1,0.25)), v=((seq(0, maxtime, myby))),
           col="lightgray", lty="dotted", lwd = 2,
           usr= c(0,maxtime,0.4,1) )
}
myplot <- function(mysurv, maxtime = tail(mysurv$time,1), myby = findmyby(maxtime)) {
    par(new=T)
    if(!exists('genelabel'))
        genelabel <- paste(gene, collapse = " and ")
    if(!exists('tunit'))
        tunit <- 'days'
    if(!exists('cancer', where = .GlobalEnv, mode = "character"))
        cancer <- mycancer

    plot(mysurv, col = mycol ,axes=F, ylim=c(0,1.05), yaxs='i', xlim=c(0, ceiling(maxtime*1.05)))
    par(tcl= -0.2)
    axis(1, at=seq(0, ceiling(maxtime*1.05), myby/5),labels=F, lwd=1, lwd.ticks=1)
    axis(2, at=seq(0, 1,     0.25/5), labels=F, lwd=1, lwd.ticks=1)
    par(tcl= -0.5)
    axis(1, at=seq(0, ceiling(maxtime*1.05), myby), lwd=0, lwd.ticks=2)
    axis(2, at=seq(0, 1,           0.25), lwd=0, lwd.ticks=2)

    title(paste0("Cancer: ", cancer, "\nGene: ", genelabel),
          xlab = paste0("Time (", tunit,")"),
          ylab = "Survival probability")
}

addlegend <- function(pos = "bottomleft", maxtime = max(mysurv$time)) {
    if(!exists('genelabel'))
        genelabel <- paste(gene, collapse = " and ")
    if( length(gene) == 2 ) {
        legend(pos, inset=c(.02, .02),box.lwd = 0,box.col = "white",
               c(paste0("High ", gene[1], " and high ", gene[2], ", n = ", nhighhigh),
                 paste0("High ", gene[1], " and low  ",  gene[2], ", n = ", nhighlow),
                 paste0("Low  ",  gene[1], " and high ", gene[2], ", n = ", nlowhigh),
                 paste0("Low  ",  gene[1], " and low  ",  gene[2], ", n = ", nlowlow)),
               fill=mycol, border = mycol, bg = "white" )
    } else {
        legend(pos, inset=c(.02, .02),box.lwd = 0,box.col = "white",
               c(paste0("High ", genelabel, ", n = ", nhigh),
                 paste0("Low  ", genelabel, ", n = ", nlow )),
               fill=mycol, border = mycol, bg = "white" )
    }
    if(pos == "bottomleft") {
        text(paste0("p-value = ", pval), x= ceiling(maxtime*1.05), y= 1, adj = c(1,1))
    } else if(pos == "topright") {
        text(paste0("p-value = ", pval), x= ceiling(maxtime*1.05) * .02, y= .05, adj = c(0,0))
    }
}

isFALSE <- function(variable) { identical(variable, FALSE) }

findtimeunit <- function(maxtime = tail(surv$Time, 1), forceunit = FALSE) {
    if (forceunit == 'years' || (isFALSE(forceunit) && maxtime > 4000)) {
        surv$Time <<- surv$Time / 365
        maxtime <<- tail(surv$Time, 1)
        return('years')
    } else if (forceunit == 'months' || (isFALSE(forceunit) && maxtime > 365)) {
        surv$Time <<- surv$Time / 30.5
        maxtime <<- tail(surv$Time, 1)
        return('months')
    } else {
        return('days')
    }
}

checkforemptygroups <- function (..., cnc) {
    groups <- c(...)
    if(sum(sign(groups)) < 2) {
        errtext <- paste0("Cancer: ", cnc, "\nGene: ", paste(gene, collapse = ' and '), "\nCohort could not be separated by gene expression")
        frame(); text(0.5,0.5,errtext )
        print(gsub(x=errtext,pattern = "\n", replacement = ". "))
        # stop("No patients match the selected expression pattern")
        return(FALSE)
    } else {
        return(TRUE)
    }
}

splitdata <- function(eset, gene, id) {
    group <- rep(0, length(eset$TCGA_Barcode) )
    if(length(gene) == 1) {
        exp.gene <- eset[[ gene ]]
        med <- median(exp.gene)
        sd  <- sd(exp.gene)

        low  <- name[exp.gene <= med]
        high <- name[exp.gene > med]

        group <- rep(0, length(id) )
        group[id %in%  low]  <- 'LOW'
        group[id %in% high] <- 'HIGH'

        return(group)

    } else if(length(gene) == 2) {     # Initiate cross-testing

        exp.gene1 <- eset[[ gene[1] ]]
        exp.gene2 <- eset[[ gene[2] ]]
        med1 <- median(exp.gene1)
        med2 <- median(exp.gene2)
        sd1  <- sd(exp.gene1)
        sd2  <- sd(exp.gene2)

        lowlow   <- name[exp.gene1 <= med1 & exp.gene2 <= med2]
        lowhigh  <- name[exp.gene1 <= med1 & exp.gene2 >  med2]
        highlow  <- name[exp.gene1 >  med1 & exp.gene2 <= med2]
        highhigh <- name[exp.gene1 >  med1 & exp.gene2 >  med2]

        group <- rep(0, length(surv$Barcode) )
        group[id %in%   lowlow] <-   'LOWLOW'
        group[id %in%  lowhigh] <-  'LOWHIGH'
        group[id %in%  highlow] <-  'HIGHLOW'
        group[id %in% highhigh] <- 'HIGHHIGH'

        return(group)

    } else {
        stop("Cross-testing only supports 2 paired genes.")
    }
}

confounders <- function(clin, geneID) {
    vars <- head(names(clin)[-1], -1)  # Remove 'barcode' and 'group'
    confs <- character(0)
    genes <- unlist(strsplit(geneID, ',', fixed = TRUE))

    if(length(genes) == 1) {
        ## Confounder-test for single gene
        lowclin  <- subset(clin, Group == 'LOW')
        highclin <- subset(clin, Group == 'HIGH')

        for(var in vars) {
            l <- table( lowclin[[var]])
            h <- table(highclin[[var]])

            count <- length(lowclin[[var]]) + length(highclin[[var]])
            count.noNA <- length(na.omit(lowclin[[var]])) + length(na.omit(highclin[[var]]))
            if (count.noNA < count/4) {
                #                confs <- append(confs, paste(var, "was skipped due to NAs"))
                next()
            }

            if ( length(l) < 2 || length(h) < 2 ) {
                next()
            }
            if ( !is.factor(clin[[var]]) ) {
                ### If the category is continuous
                if( sum(!is.na(highclin[[var]])) > 1 || sum(!is.na(lowclin[[var]])) > 1 ) {
                    welch <- t.test(lowclin[[var]], highclin[[var]])
                    pval <- signif(welch$p.value, digits = 6)

                    if(pval <= 0.05) {
                        plot(lowclin[[var]], highclin[[var]])
                        browser()
                        confs <- append(confs, paste0(var, ", pval: ", pval))
                    }
                } else {
                    browser()
                    confs <- append(confs, paste(var, "was skipped due to insufficient data"))
                }
            } else {
                ### If the category is categorical, perform fischer's exact test.
                if( !identical(names(l),names(h)) ) {
                    cat("Wrong contable\n")
                    print(l)
                    print(h)
                    confs <- append(confs, paste("ERROR.",var,"was skipped due to wrong data"))
                    next()
                }

                contable <- matrix(c(l, h), ncol = 2,
                                   dimnames = list(names(l), c('LOW','HIGH')) )

                if(any( dim(contable) > 15)) {
                    confs <- append(confs, paste("ERROR.",var,": data problem"))
                    next()
                }

                err <- try(fisher <- fisher.test(contable), silent = TRUE)
                if(is(err, 'try-error')) {
                    fisher <- fisher.test(contable, simulate.p.value = TRUE, B=1e5)
                }
                pval <- signif(fisher$p.value, digits = 6)
                if(pval <= 0.05) {
                    # barplot(contable, beside=TRUE, ylim=c(0,40), main=var, legend.text = TRUE)
                    confs <- append(confs, paste(var, ", pval: ", pval, sep=""))
                }
            }
        }
    } else {
        #### Confounder-test for pairwise gene expression
        lowlowclin   <- subset(clin, Group == 'LOWLOW')
        lowhighclin  <- subset(clin, Group == 'LOWHIGH')
        highlowclin  <- subset(clin, Group == 'HIGHLOW')
        highhighclin <- subset(clin, Group == 'HIGHHIGH')

        for(var in vars) {
            # cat(cancer,geneID,var,"\n",sep=" ")

            ### If too many NAs
            count <- length(lowlowclin[[var]])  +
                length(lowhighclin[[var]]) +
                length(highlowclin[[var]]) +
                length(highhighclin[[var]])
            count.noNA <- length(na.omit(lowlowclin[[var]]))  +
                length(na.omit(lowhighclin[[var]])) +
                length(na.omit(highlowclin[[var]])) +
                length(na.omit(highhighclin[[var]]))
            if (count.noNA < count/4) {
                next()
            }

            if( !is.factor(clin[[var]])) {
                ### If the category is continuous
                group<- c(rep('LOWLOW',  length(  lowlowclin[[var]])),
                          rep('LOWHIGH', length( lowhighclin[[var]])),
                          rep('HIGHLOW', length( highlowclin[[var]])),
                          rep('HIGHHIGH',length(highhighclin[[var]])))
                value<- c(  lowlowclin[[var]],  lowhighclin[[var]],
                            highlowclin[[var]], highhighclin[[var]])

                anov <- anova(lm(value ~ group))
                pval <- signif(anov$`Pr(>F)`[1], digits = 6)

                if( !is.na(pval) && pval <= 0.05) {
                    boxplot(value ~ group,col=2:5, main=paste(cancer,geneID,var,sep=", ") )
                    print( paste(geneID, ':', var, ", pval: ", pval, sep="") )
                    phoc <- pairwise.t.test(value, group, p.adj = "bonf")
                    print(phoc$p.value)
                    confs <- append(confs, paste0(var, ", pval: ", pval))
                }
            } else {
                ### If the category is categorical, perform chisq test.
                ll <- table(  lowlowclin[[var]])
                lh <- table( lowhighclin[[var]])
                hl <- table( highlowclin[[var]])
                hh <- table(highhighclin[[var]])
                contable <- matrix(c(ll, lh, hl, hh), ncol = 4,
                                   dimnames = list(names(ll), c('LOWLOW', 'LOWHIGH', 'HIGHLOW', 'HIGHHIGH')) )

                chisq <- tryCatch(chisq.test(contable),error=function(e) e, warning=function(w) invisible(w))

                if(is(chisq, 'warning')) {
                    err <- try(fisher <- fisher.test(contable, workspace = 2e4), silent = TRUE)
                    if(is(err, 'try-error')) {
                        fisher <- fisher.test(contable, simulate.p.value = TRUE, B=1e5)
                    }
                    pval <- fisher$p.value
                } else {
                    pval <- signif(chisq$p.value, digits = 6)
                }

                if( is.na(pval)) {
                    print( paste0(geneID, ':', var, ", fisher failed") )
                    browser()
                    next()
                } else if(pval <= 0.05) {
                    print(contable)
                    print( paste0(geneID, ':', var, ", pval: ", pval) )
                    barplot(contable, main = paste(cancer,geneID,var,sep=", "), beside = T)
                    confs <- append(confs, paste0(var, ", pval: ", pval))
                }
            }
        }

    }
    return(confs)
}

removeoutliers <- function(outlierfile = 'Outliers.txt', data) {
    outliers <- read.table(outlierfile, sep="\t", header=FALSE, colClasses='character')
    outliers = as.vector(as.matrix(outliers))

    index <- match(x = outliers, table = data[[1]], nomatch = FALSE)
    if(any(index)) {
        return(data[-index])
    } else {
        return(data)
    }
}

## OUTLIER REMOVAL: LARS METODE
removeoutliers2 <- function(outlierfile = 'Outliers.txt', data) {
    outliers <- read.table(outlierfile, sep="\t", header=FALSE, colClasses='character')
    outliers = as.vector(as.matrix(outliers))

    data <- data[!rownames(data) %in% outliers,]
    return(data)
}
