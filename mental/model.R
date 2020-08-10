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
#             phenotype - file containing types of phenotype variables for selection of regression model
#             covardata - covariate data files including one or more covariate columns (default values: NULL)
#             outfile - output file of analysis
#
# Example:    Rscript model.R testdata/test 1:1000 testdata/test_25000_pheno.txt pheno_1 testdata/test_25000_pheno_types.txt testdata/test_25000_covar.txt testdata/test_1-1000_pheno_1.txt

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
library(MASS)
library(nnet)
source("readplink.R")
options(stringsAsFactors = FALSE)

pheno <- read.table(phenodata, header=T, sep=",", strip.white=T, as.is=T)
colnames(pheno)[1] <- 'IID'
pheno <- pheno[,c('IID',sub('^','X',sub('-','.',phenoname)))]
colnames(pheno) <- c('IID','pheno')
pheno <- pheno[!is.na(pheno$pheno),]
pheno <- pheno[pheno$pheno >= 0,]

phenomap <- read.table(phenotype, header=T, strip.white=T, as.is=T)
vartype = phenomap$variable_type[phenomap$field_id==unlist(strsplit(phenoname, "-"))[1]][1]

#covar <- read.table(covardata, header=T, sep=',', strip.white=T, as.is=T)
#colnames(covar)[1] <- 'IID'
#covar <- na.omit(covar)
#covar[,3] <- factor(covar[,3])
#covar <- covar[,-3]
#covar <- covar[,-2]

if (grepl(":", snplist, fixed = TRUE)==TRUE) {
    from = as.integer(unlist(strsplit(snplist, ":"))[1])
    to = as.integer(unlist(strsplit(snplist, ":"))[2])
    snplist = from:to
} else if (grepl(",", snplist, fixed = TRUE)==TRUE) {
    snplist = as.integer(unlist(strsplit(snplist, ",")))
} else if (grepl(";", snplist, fixed = TRUE)==TRUE) {
    snplist = as.integer(unlist(strsplit(snplist, ";")))
}

geno_snps <- read.table(paste0(genodata,'.bim'), header=F, strip.white=T, as.is=T)
geno_inds <- read.table(paste0(genodata,'.fam'), header=F, strip.white=T, as.is=T)
nsamples = nrow(geno_inds)
geno_snpStats <- get_bed_geno(genodata, snplist, nsamples)

if (vartype == 'binary') {
    write(paste('SNP','CHR','BP','PVAL','A1','A2','FRQ','Z','OR','SE','N','NCASE','NCONTROL', sep='\t'), file=outfile)
}else {
    write(paste('SNP','CHR','BP','PVAL','A1','A2','FRQ','Z','BETA','SE','N', sep='\t'), file=outfile)
}
for (i in 1:ncol(geno_snpStats)) {
    snpname = geno_snps$V2[snplist][i]
    chr = geno_snps$V1[geno_snps$V2==snpname][1]
    bp = geno_snps$V4[geno_snps$V2==snpname][1]
    a1 = geno_snps$V5[geno_snps$V2==snpname][1]
    a2 = geno_snps$V6[geno_snps$V2==snpname][1]

    geno <- data.frame(geno_inds[,2], geno_snpStats[,i])
    colnames(geno) <- c('IID', 'geno')
    geno$geno <- as.integer(geno$geno)
    geno <- geno[geno$geno!=-1,]

    # merge geno and pheno
    dat <- merge(geno, pheno, by="IID", all.x=F, all.y=F, sort=F)
    # merge dat and covar
    #dat <- merge(dat, covar, by="IID", all.x=F, all.y=F)
    # remove rows with na
    dat <- na.omit(dat)
    rownames(dat) <- dat$Row.names
    dat <- dat[,-1]

    n_ind = nrow(dat)

    maf = sum(dat$geno)/(n_ind*2)
    if (maf > 1-maf) {
        a3 = a1
        a1 = a2
        a2 = a3
        dat$geno <- 2 - dat$geno
        maf = 1-maf
    }
    maf = round(maf,6)

    if (vartype == 'continuous') {
        # linear regression
        fmodel <- glm(pheno ~ ., data=dat)
        nmodel <- glm(pheno ~ 1, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[2,1]
        se = summ$coefficients[2,2]
    }else if (vartype == 'binary') {
        # logistic regression
        dat$pheno <- factor(dat$pheno)
        fmodel <- glm(pheno ~ ., family=binomial, data=dat)
        nmodel <- glm(pheno ~ 1, family=binomial, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[2,1]
        se = summ$coefficients[2,2]
        or = exp(beta)
        or = round(or,6)
        se = or * se
        n_ctrl = length(dat$pheno[dat$pheno %in% c('0')])
        n_case = length(dat$pheno[dat$pheno %in% c('1')])
    }else if (vartype == 'categorical' || vartype == 'multiple_response') {
        # multinomial logistic regression
        dat$pheno <- factor(dat$pheno)
        invisible(capture.output(fmodel <- multinom(pheno ~ ., data=dat)))
        invisible(capture.output(nmodel <- multinom(pheno ~ 1, data=dat)))
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[1,2]
        se = summ$standard.errors[1,2]
    }else if (vartype == 'ordinal') {
        # ordinal logistic regression
        dat$pheno <- factor(dat$pheno)
        fmodel <- polr(pheno ~ ., data=dat)
        nmodel <- polr(pheno ~ 1, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(Chi)"
        suppressMessages(summ <- summary(fmodel))
        beta = summ$coefficients[1,1]
        se = summ$coefficients[1,2]
    }else if (vartype == 'count') {
        # poisson regression
        fmodel <- glm(pheno ~ ., family=poisson, data=dat)
        nmodel <- glm(pheno ~ 1, family=poisson, data=dat)
        anov = anova(nmodel, fmodel, test = 'Chisq')
        pval = anov$"Pr(>Chi)"
        summ = summary(fmodel)
        beta = summ$coefficients[2,1]
        se = summ$coefficients[2,2]
    }
    if (vartype %in% c('continuous','binary','categorical','multiple_response','ordinal','count')) {
        pval = as.numeric(unlist(strsplit(as.character(pval), " "))[2])
        if (pval < 1e-99) {
           pval = 1e-99
        }
        z = qnorm(pval/2)
        if (beta < 0) {
           z = -z
        }
        pval = formatC(pval, format = "e", digits = 2)
        z = round(z,6)
        beta = round(beta,6)
        se = round(se,6)
        if (vartype == 'binary') {
            write(paste(snpname, chr, bp, pval, a1, a2, maf, z, or, se, n_ind, n_case, n_ctrl, sep='\t'), file=outfile, append=TRUE)
        }else {
            write(paste(snpname, chr, bp, pval, a1, a2, maf, z, beta, se, n_ind, sep='\t'), file=outfile, append=TRUE)
        }
    }else {
        write("pheno type does not exist", file=outfile, append=TRUE)
    }
}
write('Done', file=outfile, append=TRUE)
