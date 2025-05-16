.libPaths( c('.Rpackages',.libPaths() ) )

#'# Regional DNA methylation analysis using DMRcate
#' Using data preprocessed in our script:  
#' meth01_process_data.R

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(stringi))
suppressMessages(library(limma))
suppressMessages(library(minfi))
suppressMessages(library(missMethyl))
suppressMessages(library(ChAMP))
suppressMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
options(warn=0)

#' load the data
load("./data/processed.rda")

#'# Introduction to differential variability analysis
#' see [Phipson and Oshlack. Genome Biol 2014](https://pubmed.ncbi.nlm.nih.gov/25245051/). 
#' Impute missing Beta-values (varFit will produce an error with missingness)
sum(is.na(beta))

# the champ.impute function prints a lot of information when running, so we'll use this function to suppress the progress messages
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}
betas.impute <- hush(champ.impute(beta=beta, pd=pheno, k=5, ProbeCutoff=0.2, SampleCutoff=0.1))
betas.impute <- betas.impute$beta

#' coef parameter in varFit states which columns of design matrix correspond to the intercept and variable of interest
#' If Beta-values are used, a lofit transformation is performed within the varFit function
model <- model.matrix( ~smoker+Sex+CD4+CD8+NK+B+MO,data=pheno)
head(model)
EWAS.diffVar <- varFit(betas.impute, design=model, coef=c(1,2))
Top.diffVar <- topVar(EWAS.diffVar, coef=2, number=10, sort =TRUE)
Top.diffVar

#' Bind results with annotation
Annot <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
Annot <- subset(Annot, Annot$Name %in% manifest$ilmn_id)
Annot$ID <- substr(Annot$Name, 1, 10)
Top.diffVar$ID <- row.names(Top.diffVar)
Top.diffVar <- dplyr::left_join(Top.diffVar, Annot, by = "ID")

#' Order by chr and chromosomal position
Top.diffVar$chr <- as.numeric(gsub("chr", "", Top.diffVar$chr))
Top.diffVar[,c(2,6)] <- round(Top.diffVar[,c(2,6)],2)
Top.diffVar[order(Top.diffVar$chr, Top.diffVar$pos),c(2,5,7,8,9)]

#' We can see differences in variability among the top sites
cpg1 <- data.frame(id = colnames(betas.impute), beta = betas.impute[rownames(betas.impute) == 'cg03459282',])
cpg1 <- merge(cpg1, pheno[,c('gsm', 'smoker')], by.x = 'id', by.y = 'gsm') 
cpg1$smoker <- as.numeric(factor(cpg1$smoker))

cpg2 <- data.frame(id = colnames(betas.impute), beta = betas.impute[rownames(betas.impute) == 'cg04300100',])
cpg2 <- merge(cpg2, pheno[,c('gsm', 'smoker')], by.x = 'id', by.y = 'gsm') 
cpg2$smoker <- as.numeric(factor(cpg2$smoker))

boxplot(cpg1$beta ~ cpg1$smoker, col = c("blue", "red"), outline = F, xlab = 'smoking', ylab = 'Beta-value', names = c("non-smoker", "smoker"));points(jitter(cpg1$smoker, amount = 0.1), cpg1$beta, pch = 16)

boxplot(cpg2$beta ~ cpg2$smoker, col = c("blue", "red"), outline = F, xlab = 'smoking', ylab = 'Beta-value', names = c("non-smoker", "smoker"));points(jitter(cpg2$smoker, amount = 0.1), cpg2$beta, pch = 16)

#' Note that diffVar results may be influenced by outliers, and therefore it is helpful to visually inspect results.

rm(Top.diffVar,EWAS.diffVar,betas.impute,cpg1,cpg2);gc()


#'# Introduction to regional analyses with DMRcate
#' see [Peters et al. Bioinformatics 2015](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-8-6).  
suppressMessages(library(DMRcate)) # Popular package for regional DNA methylation analysis

#' Let's run the regional analysis using the Beta-values from our preprocessed data
#' First reorder matrix and revert to EPICv2 probe names
beta <- beta[manifest$probe_id, ]
row.names(beta) <- manifest$ilmn_id

myannotation <- cpg.annotate("array", na.omit(beta), analysis.type="differential",arraytype="EPICv2",
                             epicv2Filter = "mean", 
                             what="Beta",design=model, coef=2)
# set analysis.type = "diffVar" to return differentially variable methylated regions     
# DMRcate needs complete data. We can use na.omit() to remove rows with missing data or use the betas.impute object
# DMRcate uses the a different annotation and remaps probes with possible off-targets; therefore, some position replicates are merged                     

#' Regions are now agglomerated from groups of significant probes 
#' where the distance to the next consecutive probe is less than lambda nucleotides away
dmrcoutput.smoking <- dmrcate(myannotation, lambda=1000, C=2)

#' Let's look at the results
results.ranges <- extractRanges(dmrcoutput.smoking, genome = "hg38")
results.ranges

#' Plot the DMR using the Gviz

#' If you are interested in plotting genomic data the Gviz is extremely useful
#' Let's look at the fourth region
results.ranges[2]

# set up the grouping variables and colours
pheno$smoker<-as.factor(pheno$smoker)
cols = c("magenta","red")[pheno$smoker]
names(cols) = levels(pheno$smoker)[pheno$smoker]

#' Draw the plot a DMR\
#+ fig.width=8, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=4, CpGs=beta, phen.col=cols, what = "Beta",
         arraytype = "EPICv2", genome="hg38")

#' Clean data
rm(dmrcoutput.smoking,myannotation);gc()


#'# Introduction to regional analyses with comb-p
#' see [Pedersen et al. Bioinformatics 2012](https://academic.oup.com/bioinformatics/article/28/22/2986/240603).  
suppressMessages(library(ENmix)) 

#' Input is a data frame with colnames "chr", "start", "end", "p" and "probe"
#' Using p-values from limma
EWAS.limma <- eBayes(lmFit(beta, design=model))
EWAS.limma <- topTable(EWAS.limma, coef=2, number=Inf, sort.by="p")
EWAS.limma <- EWAS.limma[!is.na(EWAS.limma$logFC),]
#' Add chr and pos
EWAS.limma <- cbind(EWAS.limma, chr = Annot[,"chr"][match(rownames(EWAS.limma), rownames(Annot))])
EWAS.limma <- cbind(EWAS.limma, pos = Annot[,"pos"][match(rownames(EWAS.limma), rownames(Annot))])
EWAS.limma$end <- EWAS.limma$pos
#' keep only needed columns
EWAS.combp <- EWAS.limma[,c("chr", "pos", "end", "P.Value")]
colnames(EWAS.combp) <- c("chr", "start", "end", "p")
EWAS.combp$probe = rownames(EWAS.combp)
EWAS.limma$chr <- as.numeric(gsub("chr", "", EWAS.limma$chr)) # chr as numeric
# Run combp. Results saved in working directory. For efficiency, we'll use output already generated.
# combp(EWAS.combp, region_plot = F, mht_plot = F, dist.cutoff=1000, seed=0.001, verbose = T) 
# dist.cutoff = maximum distrance between basepairs to combine adjacent DMRs; seed = FDR significance threshold for initial selection of DMRs

#' read output
combp.result = read.csv('./results/resu_combp.csv')
combp.result[1:5, ]
# Regions chr1: 92946131-92947588, chr2: 233284661-233285454, chr2: 27665079-27665711, chr3: 21792434-21793157, chr5: 373299-373887, chr6: 32120773-32121261, and chr14: 106329158-106330538 overlap with DMRcate results. 
# comb-p has higher power to detect small effect sizes, but increased Type I error; see [Malid et al. Briefings in Bioinformatics 2019](https://academic.oup.com/bib/article/20/6/2224/5096828)

#' Other popular options for conducting Regional DNA methylation analysis in R are Aclust and bumphunter

#' Visualizing DMPs and DMRs using EWASplot
suppressMessages(library(ggplot2))
suppressMessages(library(EWASplot))
manhattan_plot(probe = EWAS.limma, region = combp.result)


#' clear environment
rm(list = ls()); gc()

