#--------------------------- Description ----------------------------#

# Function: This script 

# Authors: Yunhan Chu, Alexey A. Shadrin

# (c) 2020-2022 NORMENT, UiO

# Usage:      Rscript test.R genodata phenodata phenoname
#
# Arguments:  genodata  - prefix of genotype plink files
#             snplist   - a numeric or a character vector indicating a subset of SNPs to be selected (default value: NULL)
#             phenodata - phenotype data files including one or more phenotype columns
#             phenoname - a column name of phenotype in phenodata
#             phenotype - file containing types of phenotype for selection of regression model
#             covardata - covariate data files including one or more covariate columns (default values: NULL)
#             outfile - output file of analysis
#
# Example:    Rscript model.R testdata/test 1:1000 testdata/test_25000_pheno.txt pheno_1 testdata/test_25000_pheno_types.txt testdata/test_25000_covar.txt testdata/test_pheno_1_mostest.txt

#-------------------------- Input paramters -------------------------#

# Import arguments from command line
args <- commandArgs(TRUE)

if (length(args) < 6) {
  stop("Six arguments must be supplied")
}

# genotype data
genodata = args[1]

# snp list
snplist = args[2]

# phenotype data
phenodata = args[3]

# phenotype name
phenoname = args[4]

# phenotype file
phenotype = args[5]

# covariate data
covardata = args[6]

# covariate data
outfile = args[7]

#----------------------------- Start code ---------------------------#
suppressMessages(library(snpStats))
library(MASS)
library(nnet)

pheno <- read.table(phenodata, header=T, strip.white=T, as.is=T)
pheno <- pheno[,c('IID',phenoname)]
colnames(pheno) <- c('IID','pheno')

phenomap <- read.table(phenotype, header=T, strip.white=T, as.is=T)
phenotype = phenomap$TYPE[phenomap$PHENO==phenoname][1]

covar <- read.table(covardata, header=T, strip.white=T, as.is=T)
covar$Sex[covar$Sex=='Female'] <- 2
covar$Sex[covar$Sex=='Male'] <- 1
covar$Sex[!covar$Sex %in% c(1,2)] <- NA
covar <- covar[,2:8]

if (grepl(":", snplist, fixed = TRUE)==TRUE) {
    from = as.integer(unlist(strsplit(snplist, ":"))[1])
    to = as.integer(unlist(strsplit(snplist, ":"))[2])
    snplist = from:to
} else if (grepl(",", snplist, fixed = TRUE)==TRUE) {
    if (grepl("rs", snplist, fixed = TRUE)==TRUE) {
        snplist = unlist(strsplit(snplist, ","))
    } else {
        snplist = as.integer(unlist(strsplit(snplist, ",")))
    }
} else if (grepl(";", snplist, fixed = TRUE)==TRUE) {
    if (grepl("rs", snplist, fixed = TRUE)==TRUE) {
        snplist = unlist(strsplit(snplist, ";"))
    } else {
        snplist = as.integer(unlist(strsplit(snplist, ";")))
    }
} else {
    if (grepl("rs", snplist, fixed = TRUE)==TRUE) {
        snplist = unlist(strsplit(snplist, " "))
    } else {
        snplist = as.integer(unlist(strsplit(snplist, " ")))
    }
}

data_snpStats <- read.plink(genodata, select.snps = snplist)
geno_snpStats <- as(t(data_snpStats$genotypes), 'numeric')

write(paste('SNP', 'PVAL', sep='\t'), file=outfile)
for (i in 1:nrow(geno_snpStats)) {
    snpname = rownames(geno_snpStats)[i]
    geno <- data.frame(colnames(geno_snpStats), sub(".* ","", geno_snpStats[i,]))
    colnames(geno) <- c('IID', 'geno')

    # merge geno and pheno
    dat <- merge(geno, pheno, by= "IID", all.x=F, all.y=F, sort=F)
    # merge dat and covar
    dat <- merge(dat, covar, by= "IID", all.x= F, all.y= F)
    # remove rows with na
    dat <- na.omit(dat)
    rownames(dat) <- dat$Row.names
    dat <- dat[,-1]

    if (phenotype == 'continuous') {
        # linear regression
        fmodel <- glm(pheno ~ ., data=dat)
        nmodel <- glm(pheno ~ 1, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
    }else if (phenotype == 'binary') {
        # logistic regression
        dat$pheno[dat$pheno=='A'] <- 0
        dat$pheno[dat$pheno=='B'] <- 1
        dat$pheno <- as.integer(dat$pheno)
        fmodel <- glm(pheno ~ ., family=binomial, data=dat)
        nmodel <- glm(pheno ~ 1, family=binomial, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
    }else if (phenotype == 'multinomial') {
        # multinomial logistic regression
        invisible(capture.output(fmodel <- multinom(pheno ~ ., data=dat)))
        invisible(capture.output(nmodel <- multinom(pheno ~ 1, data=dat)))
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(Chi)"
    }else if (phenotype == 'ordered') {
        # ordinal logistic regression
        dat$pheno <- factor(dat$pheno)
        fmodel <- polr(pheno ~ ., data=dat)
        nmodel <- polr(pheno ~ 1, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(Chi)"
    }else if (phenotype == 'count') {
        # poisson regression
        fmodel <- glm(pheno ~ ., family=poisson, data=dat)
        nmodel <- glm(pheno ~ 1, family=poisson, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
    }
    pval = formatC(as.numeric(unlist(strsplit(as.character(pval), " "))[2]),3, format = "e", digits = 2)
    print(paste(snpname, pval, sep='    '))
    write(paste(snpname, pval, sep='\t'), file=outfile, append=TRUE)
}
