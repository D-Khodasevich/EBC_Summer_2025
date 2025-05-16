.libPaths( c('.Rpackages',.libPaths() ) )

#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  

# load the data
load("data/processed.rda")

#' load packages
options(warn=-1)
suppressPackageStartupMessages({
  library(DMRcate) # for regional analysis
  library(magrittr)
  library(data.table) # for fast aggregation of large data 
  library(qqman) # for visualization of data
  library(stringi) # string manipulation
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})
options(warn=0)

#' Code categorical variable as factors
pheno$smoker %<>% factor
pheno$Sex    %<>% factor


#'# Running an Epigenome Wide Association
#' Here we run an EWAS on Smoking status (as a predictor of methylation)  
nCpG = dim(beta)[1]
nCpG

#' First we can run a linear regression on a single CpG that we have already picked
CpG.name = "cg05575921"
pheno$CpG.level <- beta[CpG.name,]

#' Difference in methylation between smokers and non-smokers for this CpG
#' some descriptive statistics

pheno[,.(
   Min    = min   (CpG.level) %>% round(3)
  ,Mean   = mean  (CpG.level) %>% round(3)
  ,Median = median(CpG.level) %>% round(3)
  ,Max    = max   (CpG.level) %>% round(3)
  ,SD     = sd    (CpG.level) %>% round(3)
  ,.N),by=smoker] %>% knitr::kable(.)

#' Difference in beta methylation values between Smokers and non smokers
boxplot(CpG.level ~ smoker,pheno,main=paste0("Beta-values\n", CpG.name), col=c("blue","red"),ylim=c(.3,1))

#' Linear regression on betas
lm(CpG.level ~ smoker,data=pheno) %>% summary %>% coef

#' Comparison with m-values
pheno[,CpG.mlevel:=log2(CpG.level/(1-CpG.level))]

pheno[,.(
   Min    = min   (CpG.mlevel) %>% round(3)
  ,Mean   = mean  (CpG.mlevel) %>% round(3)
  ,Median = median(CpG.mlevel) %>% round(3)
  ,Max    = max   (CpG.mlevel) %>% round(3)
  ,SD     = sd    (CpG.mlevel) %>% round(3)
  ,.N),by=smoker] %>% knitr::kable(.)


par(mfrow=c(1,2))
boxplot(CpG.level  ~ smoker,data=pheno,main=paste0("Beta-values\n",CpG.name), col=c("blue","red"))
boxplot(CpG.mlevel ~ smoker,data=pheno,main=paste0("M-values\n"   ,CpG.name), col=c("blue","red"))
dev.off()

#' linear regression on m-values
lm(CpG.mlevel ~ smoker,data=pheno) %>% summary %>% coef

#' We can always extract measures of the relative quality of statistical models - e.g. adjusted R2 - to look at model performance  
#' model on betas
lm(CpG.level  ~ smoker,data=pheno) %>% summary %$% adj.r.squared

#' model on mvalues
lm(CpG.mlevel ~ smoker,data=pheno) %>% summary %$% adj.r.squared






#' # EWAS using limma
#' see [Smyth GK. Stat Appl Genet Mol Biol 2004](https://www.ncbi.nlm.nih.gov/pubmed/16646809).  
suppressMessages(library(limma,minfi))

#' Smoking as predictor  
#' note that limma is quite fast for running almost a million regressions!
pheno[,smoke_dummy:=ifelse(smoker=="smoker",1,0)]

#' First we need to define a model
model <- model.matrix( ~smoker,data=pheno) # simple bivariate model
EWAS.limma_model1 <- eBayes(lmFit(beta, design=model))
results1 <- topTable(EWAS.limma_model1, coef=2, number=Inf, sort.by="p", adjust.method="BH")
results1$ID <- row.names(results1)

#' here effect size is logFC among smokers compared to non-smokers for the top hits
head(results1[, c(1, 2, 5)])

#' check with previous result on our selected CpG (running lm without CpGassoc)
subset(results1, results1$ID == CpG.name)
summary(lm(CpG.level~smoker,pheno))$coefficients[2, ]
#' The effect estimates are identical
#' However the test statistics and p-values differ slightly due to 
#' eBayes usage in the limma models

#' Bonferroni significant hits
table(results1$P.Value < 0.05/(nCpG))
#' FDR significant hits
table(results1$adj.P.Val < 0.05)


#' EWAS with adjustment for cell types
#' now we can run the linear regression on betas adjusting for cell proportions
#' Need sex as indicator in covariate matrix
model <- model.matrix( ~smoker+Sex+CD4+CD8+NK+B+MO,data=pheno)

EWAS.limma_model2 <- eBayes(lmFit(beta, design=model))
results2 <- topTable(EWAS.limma_model2, coef=2, number=Inf, sort.by="p", adjust.method="BH")
results2$ID <- row.names(results2)
head(results2[, c(1, 2, 5)])


#' Using mvalues
#' Use of M-values reduces heteroscedasticity to meet linear model assumptions, 
#' See [Du P, et al. BMC Bioinformatics. 2010](https://pubmed.ncbi.nlm.nih.gov/21118553/). 
mvals <- B2M(beta) # convert beta matrix to m-value matrix
EWAS.limma_model3 <- eBayes(lmFit(mvals, design=model))
results3 <- topTable(EWAS.limma_model3, coef=2, number=Inf, sort.by="p", adjust.method="BH")
head(results3[, c(1, 2, 5)])


#' ## Genomic inflation in EWAS
#' qqplot and lambda interpretation  
#+ fig.width=13, fig.height=7, dpi=300
par(mfrow=c(1,2))
qqman::qq(results1$P.Value, main = "QQ plot for association between \n methylation and Smoking \n not adjusted for cell proportions") 
qqman::qq(results2$P.Value, main = "QQ plot for association between \n methylation and Smoking \n adjusted for cell proportions") 
dev.off()

#' Lambda - this is a summary measure of genomic inflation  
#' ratio of observed vs expected median p-value - is there early departure of the qqline
#' estimated at -log10(median=0.5) ~ 0.3 on the x-axis of a qqplot  
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

#' Lambda before cell type adjustment
lambda(results1$P.Value)
#' Lambda after cell type adjustment
lambda(results2$P.Value)
#' Lambda after cell type adjustment and using m-values
lambda(results3$P.Value)


#' Bind results with annotation
Annot <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
Annot <- subset(Annot, Annot$Name %in% manifest$ilmn_id)
Annot$ID <- substr(Annot$Name, 1, 10)
results2 <- dplyr::left_join(results2, Annot, by = "ID")
results2 %>% arrange(P.Value) %>% dplyr::select(ID, logFC, P.Value, UCSC_RefGene_Name) %>% head(n = 10)


#' Volcano Plot-results2
#' Bonferroni threshold
#' 
plot(-log10(P.Value) ~ logFC, data=results2, xlab="logFC", ylab="-log10(p-value)",
     main="Volcano Plot\nadjusted for cell proportions",xlim = c(-1,1), ylim=c(0,15))
abline(h = -log10(0.05/(nCpG)), lty=1, col="#FDE725FF", lwd=2)

#'## Manhattan plot for cell-type adjusted EWAS  
results2$chr <- as.numeric(gsub("chr", "", results2$chr))
qqman::manhattan(results2,chr="chr", bp="pos", p="P.Value", snp="ID"
   ,suggestiveline=FALSE, genomewideline = -log10(0.05/(nCpG)), ylim=c(0,15)
   ,main = "Manhattan Plot \n adjusted for cell proportions")




#'# Introduction to pathway analysis with gometh
#' see [Geeleher et al. Bioinformatics 2013](https://academic.oup.com/bioinformatics/article/29/15/1851/265573?login=false).
#' Using all FDR significant associations from the cell-adjusted model as input
library(missMethyl)
keggResult = gometh(sig.cpg = results2$ID[results2$adj.P.Val < 0.05], 
                    all.cpg = results2$ID, collection = 'KEGG')
keggResult = keggResult[order(keggResult$P.DE),]
head(keggResult)




#' ## Other tools for running EWAS
#' CpGassoc
#' See [Barfield et al. Bioinformatics 2012](http://www.ncbi.nlm.nih.gov/pubmed/22451269)  
#' See previous version of the bootcamp for examples 
#' https://github.com/D-Khodasevich/EBC_Summer_2024

#' Efficient EWAS using meffil
#' See [meffil GitHub page](https://github.com/perishky/meffil)
#' Options for running robust linear regression
#' ewas.rlm<-meffil.ewas(beta, variable=pheno$smoker, covariates=pheno[,c('Sex','CD4','CD8','NK','MO','GR','B')], rlm = TRUE)

#' clear environment
rm(list = ls()); gc()
#' End of script 03
