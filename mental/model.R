#--------------------------- Description ----------------------------#

# Function: This script runs association test between genotype and phenotype data using matched regression model.

# Authors: Yunhan Chu, Alexey A. Shadrin

# (c) 2020-2022 NORMENT, UiO

# Usage:      Rscript model.R genodata snplist phenodata phenoname phenotype covardata outfile
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
options(stringsAsFactors = FALSE)

pheno <- read.table(phenodata, header=T, strip.white=T, as.is=T)
pheno <- pheno[,c('IID',phenoname)]
colnames(pheno) <- c('IID','pheno')

phenomap <- read.table(phenotype, header=T, strip.white=T, as.is=T)
phenotype = phenomap$TYPE[phenomap$PHENO==phenoname][1]

covar <- read.table(covardata, header=T, strip.white=T, as.is=T)
covar$Sex <- factor(covar$Sex)
covar$genBatch <- factor(covar$genBatch)
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

geno_bim <- read.table(paste0(genodata,'.bim'), header=F, strip.white=T, as.is=T)
data_snpStats <- read.plink(genodata, select.snps = snplist)
geno_snpStats <- as(t(data_snpStats$genotypes), 'numeric')

if (phenotype == 'binary') {
    write(paste('SNP','CHR','BP','PVAL','A1','A2','MAF','NCASE','NCONTROL','Z','OR','SE', sep='\t'), file=outfile)
}else {
    write(paste('SNP','CHR','BP','PVAL','A1','A2','MAF','N','Z','BETA','SE', sep='\t'), file=outfile)
}
for (i in 1:nrow(geno_snpStats)) {
    snpname = rownames(geno_snpStats)[i]
    chr = geno_bim$V1[geno_bim$V2==snpname][1]
    bp = geno_bim$V4[geno_bim$V2==snpname][1]
    a1 = geno_bim$V5[geno_bim$V2==snpname][1]
    a2 = geno_bim$V6[geno_bim$V2==snpname][1]

    geno <- data.frame(colnames(geno_snpStats), sub(".* ","", geno_snpStats[i,]))
    colnames(geno) <- c('IID', 'geno')
    geno$geno <- as.integer(geno$geno)

    # merge geno and pheno
    dat <- merge(geno, pheno, by= "IID", all.x=F, all.y=F, sort=F)
    # merge dat and covar
    dat <- merge(dat, covar, by= "IID", all.x= F, all.y= F)
    # remove rows with na
    dat <- na.omit(dat)
    rownames(dat) <- dat$Row.names
    dat <- dat[,-1]

    n_ind = nrow(dat)
    if (a2 > a1) {
        dat$geno <- 2 - dat$geno
    }
    maf = sum(dat$geno)/(n_ind*2)
    # let the A1 to be minor allele
    if (maf > 1-maf) {
        a3 = a1
        a1 = a2
        a2 = a3
        dat$geno <- 2 - dat$geno
        maf = 1-maf
    }
    maf = round(maf,6)

    if (phenotype == 'continuous') {
        # linear regression
        fmodel <- glm(pheno ~ ., data=dat)
        nmodel <- glm(pheno ~ 1, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[2,1]
        se = summ$coefficients[2,2]
    }else if (phenotype == 'binary') {
        # logistic regression
        dat$pheno <- factor(dat$pheno)
        fmodel <- glm(pheno ~ ., family=binomial, data=dat)
        nmodel <- glm(pheno ~ 1, family=binomial, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[2,1]
        se = summ$coefficients[2,2]

        n_ctrl = length(dat$pheno[dat$pheno %in% c('0','A')])
        n_case = length(dat$pheno[dat$pheno %in% c('1','B')])
        n_ctrl_a1 = sum(dat$geno[dat$pheno %in% c('0','A')])
        n_case_a1 = sum(dat$geno[dat$pheno %in% c('1','B')])
        n_ctrl_a2 = n_ctrl * 2 - n_ctrl_a1
        n_case_a2 = n_case * 2 - n_case_a1

        or = (n_case_a1 * n_ctrl_a2)/(n_ctrl_a1 * n_case_a2)
        or = round(or,6)
        se = sqrt(1/n_case_a1 + 1/n_ctrl_a1 + 1/n_case_a2 + 1/n_ctrl_a2)
    }else if (phenotype == 'multinomial') {
        # multinomial logistic regression
        dat$pheno <- factor(dat$pheno)
        invisible(capture.output(fmodel <- multinom(pheno ~ ., data=dat)))
        invisible(capture.output(nmodel <- multinom(pheno ~ 1, data=dat)))
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[1,2]
        se = summ$standard.errors[1,2]
    }else if (phenotype == 'ordered') {
        # ordinal logistic regression
        dat$pheno <- factor(dat$pheno)
        fmodel <- polr(pheno ~ ., data=dat)
        nmodel <- polr(pheno ~ 1, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(Chi)"
        suppressMessages(summ <- summary(fmodel))
        beta = summ$coefficients[1,1]
        se = summ$coefficients[1,2]
    }else if (phenotype == 'count') {
        # poisson regression
        fmodel <- glm(pheno ~ ., family=poisson, data=dat)
        nmodel <- glm(pheno ~ 1, family=poisson, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[2,1]
        se = summ$coefficients[2,2]
    }
    if (phenotype %in% c('continuous','binary','multinomial','ordered','count')) {
        pval = as.numeric(unlist(strsplit(as.character(pval), " "))[2])
        z = qnorm(pval/2)
        if (beta < 0) {
           z = -z
        }
        pval = formatC(pval, format = "e", digits = 2)
        z = round(z,6)
        beta = round(beta,6)
        se = round(se,6)
        if (phenotype == 'binary') {
            write(paste(snpname, chr, bp, pval, a1, a2, maf, n_case, n_ctrl, z, or, se, sep='\t'), file=outfile, append=TRUE)
        }else {
            write(paste(snpname, chr, bp, pval, a1, a2, maf, n_ind, z, beta, se, sep='\t'), file=outfile, append=TRUE)
        }
    }else {
        write("pheno type does not exist", file=outfile, append=TRUE)
    }
}
write('Done', file=outfile, append=TRUE)
