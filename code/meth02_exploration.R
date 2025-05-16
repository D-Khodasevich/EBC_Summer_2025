.libPaths( c('.Rpackages',.libPaths() ) )

options(warn=-1)
library(data.table)
suppressMessages(library(limma))
suppressMessages(library(ENmix))
suppressMessages(library(sva))
suppressMessages(library(ChAMP))
options(warn=0)


#+ setdir01, echo = F
# knitr::opts_knit$set(root.dir = "../")
load("data/processed.rda")

#'# Exploring global DNA Methylation variability via PCs
# Let's look at the effect of sex (note that X/Y chromosomes have been removed)
plotMDS(beta, top=10000, gene.selection="common",
        pch=17,col=c("deeppink","blue")[factor(pheno$Sex)],
        dim=c(1,2),cex=1.5)
legend("right", legend=levels(factor(pheno$Sex)),bty='n',
       cex=1.5,pch=17,col=c("deeppink","blue"))
       
# Let's look at the effect of technical variables
# Row
cols <- rainbow(n = length(levels(factor(pheno$row))))
plotMDS(beta, top=10000, gene.selection="common",
        pch=17,col=cols[factor(pheno$row)],
        dim=c(1,2),cex=1.5)
legend("right", legend=levels(factor(pheno$row)),bty='n',
       cex=1.5,pch=17,col=cols)
# If clustering by technical variables is present, batch effects can be adjusted for using ComBat or using covariates for batch in downstream analyses.
# Example of ComBat
# Impute missing Beta-values (ComBat will produce an error with missingness)
# sum(is.na(beta))
# betas.impute = champ.impute(beta=beta, pd=pheno, k=5, ProbeCutoff=0.2, SampleCutoff=0.1)
# betas.impute = betas.impute$beta
# Run ComBat
# batch <- factor(pheno$Sentrix_ID)
# modcombat <- model.matrix(~1, data = pheno)
# betaCombat <- ComBat(dat = as.matrix(betas.impute), batch = batch, mod = modcombat)


#' It would be useful to look at several traits with global variability
cov<-data.frame(pheno[,c('smoker','Sex','CD4','CD8','NK','MO','GR','B')])
npc <- 20 # Top 20 PCs
svd <- prcomp(t(na.omit(beta)))
screeplot(svd, npc, type = "barplot")
eigenvalue <- svd[["sdev"]]^2
prop <- (sum(eigenvalue[1:npc])/sum(eigenvalue)) * 100
cat("Top ", npc, " principal components can explain ", 
    prop, "% of data \n    variation", "\n")
screeplot(svd,npc,type="barplot")
pcrplot(na.omit(beta), cov, npc=10) # Already saved in working directory


#'# Epigenetic Age
#' We are using epigenetic age predictions from the methscore function in the ENmix R package
#' See: https://doi.org/10.1093/bioinformatics/btae302 
#' Usage of methscore is detailed below, we are loading in previously generated results
#suppressMessages(library(ENmix))
#' phenotype file requires specifically named Age and Female variables
#names(pheno)[names(pheno) == "age"] <- "Age"
#names(pheno)[names(pheno) == "gsm"] <- "SampleID"
#pheno$Female <- ifelse(pheno$Sex == "m", 0, 1)
#clocks <- methscore(beta, datPheno=pheno, fastImputation=FALSE, normalize=TRUE,
#           GrimAgeComponent=NULL, UserRef=NULL, ForceUserRef=FALSE)
#write.csv(clocks, "./results/methscore_predictions.csv", row.names = FALSE)

#' Check missingness of CpGs for epigenetic clocks
methscore_summary <- read.csv("./results/summary_methscore_CpG.csv")
methscore_summary[1:4, c(1, 7:9)]

#' Examine fit of some clocks with chronological age
clocks <- read.csv("./results/methscore_predictions.csv")

#' Horvath's panTissue clock
cor.test(clocks$Age,clocks$HorvathAge)
plot(clocks$Age,clocks$HorvathAge,pch=21,xlab="Chronological Age",
     ylab="Horvath's DNAm Age",cex=1.2, main="Epigenetic Clocks")
legend("topleft",legend=paste0("r=",round(cor(clocks$Age,clocks$HorvathAge),2)),bty="n")
abline(lm(clocks$HorvathAge~clocks$Age),col="red",lw=2)

#' PC GrimAge
cor.test(clocks$Age,clocks$PCGrimAge)
plot(clocks$Age,clocks$PCGrimAge,pch=21,xlab="Chronological Age",
     ylab="PC GrimAge",cex=1.2, main="Epigenetic Clocks")
legend("topleft",legend=paste0("r=",round(cor(clocks$Age,clocks$PCGrimAge),2)),bty="n")
abline(lm(clocks$PCGrimAge~clocks$Age),col="red",lw=2)

#' DNAmTL
cor.test(clocks$Age,clocks$DNAmTL)
plot(clocks$Age,clocks$DNAmTL,pch=21,xlab="Chronological Age",
     ylab="DNAmTL Telomere Length Estimator",cex=1.2, main="Epigenetic Clocks")
legend("topright",legend=paste0("r=",round(cor(clocks$Age,clocks$DNAmTL),2)),bty="n")
abline(lm(clocks$DNAmTL~clocks$Age),col="red",lw=2)

group <- NA
group[pheno$smoker==1] <- 1
group[pheno$smoker==0] <- 2
pairs(~ Age + HorvathAge + Horvath2 + HannumAge + 
        PCGrimAge + DNAmTL + PhenoAge, 
      col = c("red","blue")[group],
      pch = c(8, 18)[group], 
      data = clocks)



#' Age Acceleration Residuals
#' Calculate EAA measures as residuals from a regression of 
#' Epigenetic Age Measure ~ Chronological Age
#' Note the methscore function also calculates residuals for each clock measure
clocks$Horvath_resid <- residuals(lm(HorvathAge ~ Age, data = clocks))
clocks$PCGrimAge_resid <- residuals(lm(PCGrimAge ~ Age, data = clocks))
clocks$PCPACKYRS_resid <- residuals(lm(PCPACKYRS ~ Age, data = clocks))

#' Differences by Smoking status
boxplot(clocks$Horvath_resid ~ pheno$smoker, col=c("blue","red"))
wilcox.test(clocks$Horvath_resid ~ pheno$smoker)

boxplot(clocks$PCGrimAge_resid ~ pheno$smoker, col=c("blue","red"))
wilcox.test(clocks$PCGrimAge_resid ~ pheno$smoker)

#' Smoking packyears estimator component of PC GrimAge
boxplot(clocks$PCPACKYRS_resid ~ pheno$smoker, col=c("blue","red"))
wilcox.test(clocks$PCPACKYRS_resid ~ pheno$smoker)



#' Protein Methylation Risk Scores
#' methscore also features a number of epigenetic predictors of other traits
#' including methylation risk scores of the circulating proteome
#' https://pmc.ncbi.nlm.nih.gov/articles/PMC8880990/
proteins <- clocks[, 210:318]
for(i in 1:ncol(proteins)){
  proteins[, i] <- scale(proteins[, i])
  fit <- glm(proteins[, i] ~ pheno$smoker)
  sums <- summary(fit)$coefficients[2, ]
  sums[5] <- colnames(proteins)[i]
  if(i == 1){fullsums <- sums}
  else{fullsums <- rbind(sums, fullsums)}
}
fullsums <- as.data.frame(fullsums)
colnames(fullsums) <- c("Estimate", "SE", "t", "p", "Name")
fullsums$Estimate <- as.numeric(fullsums$Estimate); fullsums$p <- as.numeric(fullsums$p)
fullsums$sig <- ifelse(fullsums$p < 0.05, fullsums$Name, "")
library(ggplot2)
ggplot(fullsums) + 
  geom_point(aes(x = Estimate, y = -log10(p))) + theme_classic() + 
  ggrepel::geom_text_repel(aes(x = Estimate, y = -log10(p), label = sig), color = "steelblue") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  ggtitle("Assocations between Smoking Status and Each MRS")


#' There are a number of other resources for calculating epigenetic clocks
#' methylCIPHER R package: https://github.com/MorganLevineLab/methylCIPHER
#' -  Another comprehensive GitHub package
#' Horvath/Clock Foundation Online Calculator](http://dnamage.genetics.ucla.edu/)
#' -  Most common original method. Option for submitting IDATs now
#' methylclock R package: https://bioconductor.org/packages/release/bioc/html/methylclock.html
#' -  Another comprehensive R package
#' Bio-learn package: https://bio-learn.github.io/
#' -  Python-based resource for clock calculations
#' -  Also the primary public implementation of GrimAge/GrimAge2




#'# Predicting smoking with EpiSmokEr
#' see [Bollepalli, Sailalitha, et al. Epigenomics 2019](https://pubmed.ncbi.nlm.nih.gov/31466478/).  
suppressMessages(require(EpiSmokEr))
# Make sure rows of pheno match betas column names
rownames(pheno)<-pheno$gsm
identical(colnames(beta),rownames(pheno))

# pheno needs a column for sex,in the format of 1 and 2 representing men and women respectively
pheno$sex<-ifelse(pheno$sex=="m",1,2)
# 121 CpGs are used selected by LASSO along with Sex to get 3 categories (current, former and never smokers)
result <- epismoker(dataset=beta, samplesheet = pheno, method = "SSt")
# Let's look how well the prediction performed
table(pheno$smoker,result$PredictedSmokingStatus)


#' clear environment
rm(list = ls()); gc()

#' End of script 02