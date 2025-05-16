# preparatory analysis to identify smokers/non-smokers
library(GEOquery)
library(data.table)
library(sesame)
library(EpiSmokEr)
smoke_cpgs <- EpiSmokEr::Illig_data;  smoke_cpgs <- smoke_cpgs$cpgs
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno = as.data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
anno <- subset(anno, anno$Methyl450_Loci %in% smoke_cpgs)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246337
betas <- fread("GSE246337_betas.csv.gz")
betas <- betas[V1 %in% anno$Name, ]
betas <- as.matrix(betas)
cpglist <- intersect(anno$Name, betas[, 1])
anno <- subset(anno, anno$Name %in% cpglist)
row.names(betas) <- betas[, 1]
betas <- betas[, -1]
betas <- betas[anno$Name, ]
row.names(betas) <- anno$Methyl450_Loci
storage.mode(betas) <- "numeric"
result_SSc <- epismoker(dataset = betas, method = "SSc")
library(dplyr)
result_SSc %>% arrange(smokingScore) %>% head(n = 15) %>% select(SampleName) %>% as.vector()
result_SSc %>% arrange(desc(smokingScore)) %>% head(n = 15) %>% select(SampleName) %>% as.vector()

# lowest 15
lowest <- c("207700460019_R01C01", "207686150008_R02C01", "207805820146_R04C01",
 "207700460090_R04C01","207705770129_R02C01", "207686140115_R07C01", "207700460028_R01C01",
 "207805820127_R03C01","207700470041_R03C01", "207805820161_R08C01", "207700460019_R07C01",
 "207686140117_R04C01","207705770050_R04C01", "207700460019_R06C01", "207700470029_R07C01")

# highest 15
highest <- c("207700460028_R06C01", "207700460020_R03C01", "207700470036_R08C01",
  "207700460111_R01C01","207805820146_R01C01", "207700460035_R02C01", "207700460069_R04C01",
  "207700460028_R03C01","207700460068_R02C01", "207805820128_R04C01", "207700470024_R06C01",
  "207700460034_R06C01","207686140125_R02C01", "207686150104_R04C01", "207805820161_R04C01")

# sex-mismatched sample
# GSM7867094

# load phenotype file to link basename to gsm
pheno = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246337&targ=gsm&form=text&view=brief"
pheno = readLines(pheno)
pheno = split(pheno,cumsum(pheno %like% "^\\^SAMPLE = GSM"))
names(pheno) = map(pheno,1) %>% stri_match_first(regex="GSM\\d+")
imap(pheno,function(s,acc){
  s = strsplit(s,split=" = ",fixed=TRUE)	
  data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> pheno
pheno = rbindlist(pheno)
pheno = pheno[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]
i = pheno[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = pheno$value[i] %>% stri_split(fixed=": ")
pheno$variable[i] = map_chr(ch,1)
pheno$value   [i] = map_chr(ch,2)
rm(ch)
pheno[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
pheno[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]
pheno = dcast(pheno, gsm ~ variable)
pheno$basename <- substr(pheno$grn, 79, 97)
pheno <- subset(pheno, basename %in% highest | basename %in% lowest | gsm == "GSM7867094")
pheno %>% dplyr::filter(basename %in% highest) %>% dplyr::select(gsm) %>% as.vector()
pheno %>% dplyr::filter(basename %in% lowest) %>% dplyr::select(gsm) %>% as.vector()

# additional mismatched tissue type samples obtained from: 
# t cells
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269759
# other tissue (brain tumor)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274910


# download datasets
.libPaths( c('.Rpackages',.libPaths() ) )

library(ewastools)
library(stringi)
library(data.table)
library(magrittr)
library(purrr)

# Download phenotype data (copy URL in browser to see what the requested file looks like)
pheno = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246337&targ=gsm&form=text&view=brief"
pheno = readLines(pheno)

# Split into individual samples
pheno = split(pheno,cumsum(pheno %like% "^\\^SAMPLE = GSM"))

# Extract GSM accessions
names(pheno) = map(pheno,1) %>% stri_match_first(regex="GSM\\d+")

# Parse pheno data
imap(pheno,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> pheno

pheno = rbindlist(pheno)

# Keep only information on sample characteristics and supplementary files
pheno = pheno[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]
i = pheno[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = pheno$value[i] %>% stri_split(fixed=": ")
pheno$variable[i] = map_chr(ch,1)
pheno$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
pheno[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
pheno[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

# Reshape data.table from long to wide format
pheno = dcast(pheno, gsm ~ variable)

# Select and parse the relevant variables
pheno = pheno[,.(gsm,age,Sex,red,grn)]
pheno = rbind(pheno
	,list(
	 "GSM8326502" 
	,"66"
	,"M"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8326nnn/GSM8326502/suppl/GSM8326502_207755830088_R01C01_Red.idat.gz"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8326nnn/GSM8326502/suppl/GSM8326502_207755830088_R01C01_Grn.idat.gz")
	,list(
	 "GSM8461228"
	,"54"
	,"F"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8461nnn/GSM8461228/suppl/GSM8461228_208066090142_R07C01_Red.idat.gz"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8461nnn/GSM8461228/suppl/GSM8461228_208066090142_R07C01_Grn.idat.gz")
	 )


# Select the samples
setkey(pheno,"gsm")
pheno = pheno[c(
	## 15 smokers
  "GSM7866977", "GSM7867039", "GSM7867098", "GSM7867112", "GSM7867122", "GSM7867159", "GSM7867181",
  "GSM7867204", "GSM7867222", "GSM7867227", "GSM7867316", "GSM7867370", "GSM7867379", "GSM7867386",
  "GSM7867450",
	
	##  15 non-smokers
  "GSM7867006", 
  "GSM7867307", # random
  "GSM7867049", "GSM7867114", "GSM7867165", "GSM7867175", "GSM7867179",
  "GSM7867250", "GSM7867261", "GSM7867294", "GSM7867369", "GSM7867388", "GSM7867389", "GSM7867419",
  "GSM7867422"

  ,"GSM7867094" # sex mismatched sample  

,"GSM8326502" # unrelated sample from another GSE, T cell instead of whole blood
,"GSM8461228" # unrelated sample of brain tissue
)]

dir.create("data",showWarnings=FALSE)

# Adding chip and position information
pheno$Sentrix_ID = substring(pheno$red, 79, 90)
pheno$row = substring(pheno$red, 92, 94)
pheno$col = substring(pheno$red, 95, 97)

# Download .idat files
map2(pheno$red,"data/" %s+% pheno$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(pheno$grn,"data/" %s+% pheno$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible
pheno$red = NULL; pheno$grn = NULL

# Import the methylation data
meth = read_idats("data/" %s+% pheno$gsm)

# add predicted smoking activity
pheno$smoker <- c(rep(1, 15), rep(0, 15), 0, 0, 0)

# recode Sex
pheno$Sex = tolower(pheno$Sex)

pheno = pheno[,.(gsm,age,Sex,smoker,Sentrix_ID,row,col)]

write.csv(pheno,file="data/pheno.csv",row.names=FALSE)
write.csv(pheno[1:30],file="data/pheno_clean.csv",row.names=FALSE)
