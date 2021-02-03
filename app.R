#  hdsheep.cer.auckland.ac.nz Domain name?


#setwd("/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code")
library(shiny)
#install.packages(devtools)
library(devtools, quietly = TRUE)
#install_github('rstudio/gt')
library(exCorr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DiffCorr")
#biocLite("mixOmics")

library(mixOmics)
library(DiffCorr)
library(gridExtra)
library(grid)
library(png)
library(shinyjs)
library(DT)
library(rJava)
library(mailR)




# the cor.prob() function below creates a correlation matrix, which also returns the p-value
# i.e. with cor() a symmetrical matrix is produced with only correlation values
# with cor.prob() an unsymmetrical matrix is produced with half correlation values half p-values

cor.prob <- function(X, dfr = nrow(X) -2){
  R <- cor(X, use = "pairwise.complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <-NA
  R
}

# the flattensequareMatrix() function below takes the cor.prob() matrix, and breaks it down into 4 columns:
# rownames, column names, correlation and p-value

flattensequareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a sequare matrix.")
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor = t(m)[ut],
             p = m[ut])
}

my_compcorr <- function(x){
  compcorr(n1 = 6, r1 = x[,3], n2 = 6, r2 = x[,5])$pval
}

# Read in sheep info table for dataset page
#sheep_info = read.csv("Sheep_info.csv", stringsAsFactors = FALSE)

# Read in nanoString dataset
nano_24 <- read.csv("EM_USE THIS averaged nanoString data normalised_MG.csv", stringsAsFactors = FALSE)
# Separate full dataset into two tissue specific datasets
nano_24_DL <- nano_24[nano_24$tissue == "nanoString Dorsolateral", ]
nano_24_DL$tissue = "Dorsolateral"
# Setting row name to sample name
rownames(nano_24_DL) <- nano_24_DL$sample
# Removing sample column - unnecessary column due to row names
nano_24_DL <- nano_24_DL[-1]
# Repeat for other tissue type
nano_24_DM <- nano_24[nano_24$tissue == "nanoString Dorsomedial", ]
nano_24_DM$tissue = "Dorsomedial"
rownames(nano_24_DM) <- nano_24_DM$sample
nano_24_DM <- nano_24_DM[-1]

#############  find all T and C cors for DM tissue dataset
# Selecting measured variables from within the dataset that are used for the analysis
DM_CT_dataset <- nano_24_DM[1:12, 5:30]
# Calling two functions on the selected variables to create the dataset of all variables correlated against one another for control samples
DM_CT_control_cors <- flattensequareMatrix(cor.prob(DM_CT_dataset[1:6,])) # control cors
# Calling two functions on the selected variables to create the dataset of all variables correlated against one another for transgenic samples
DM_CT_transgenic_cors <- flattensequareMatrix(cor.prob(DM_CT_dataset[7:12,])) # transgenic cors
# Combine the two correlation datasets
DM_C_and_T_cors <- cbind(DM_CT_control_cors,DM_CT_transgenic_cors[3:4])
# Set the column names to useful names
colnames(DM_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# Calculate Fisher value for correlation statistics
DM_CT_fisher = as.array(my_compcorr(DM_C_and_T_cors))
# Append Fisher value to data frame as a 7th column
DM_C_and_T_cors_compcor = cbind(DM_C_and_T_cors, "fisher-r-to-z" = DM_CT_fisher)


############ Find all T and C cors for DL tissue dataset
# Previous steps are repeated for the other tissue type from the nanoString analysis
DL_CT_dataset <- nano_24_DL[1:12, 5:30]
DL_CT_control_cors <- flattensequareMatrix(cor.prob(DL_CT_dataset[1:6,])) # control cors
DL_CT_transgenic_cors <- flattensequareMatrix(cor.prob(DL_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
DL_C_and_T_cors <- cbind(DL_CT_control_cors,DL_CT_transgenic_cors[3:4])
colnames(DL_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column
DL_CT_fisher = as.array(my_compcorr(DL_C_and_T_cors))
DL_C_and_T_cors_compcor = cbind(DL_C_and_T_cors, "fisher-r-to-z" = DL_CT_fisher)

################ Import metabolite data #################
# Read in metabolite data and split into the four respective tissue types
metab_all <- read.csv("EM_5yrMetab_allAvg_MG.csv", stringsAsFactors = FALSE)
# Separate metabolite file by tissue type
metab_cb <- metab_all[metab_all$tissue == "GC-MS Metabolites Cerebellum", ]; rownames(metab_cb) <- metab_cb$sample
metab_hipp <- metab_all[metab_all$tissue == "GC-MS Metabolites Hippocampus", ]; rownames(metab_hipp) <- metab_hipp$sample; metab_hipp$Benzoic.acid <- NULL # Benzoic acid variable is essentially removed from the hippocampus tissue due to not being measured. The variable is kept in the dataset to ensure the datasets are the same size, but set to NA to ensure that the variable has no effect on the statistics.
metab_liv <- metab_all[metab_all$tissue == "GC-MS Metabolites Liver", ]; rownames(metab_liv) <- metab_liv$sample
metab_mctx <- metab_all[metab_all$tissue == "GC-MS Metabolites Motor Cortex", ]; rownames(metab_mctx) <- metab_mctx$sample

#############  find all T and C cors for CB tissue dataset
# The steps followed here are the same as those conducted on the nanoString datasets
CB_CT_dataset <- metab_cb[1:12, 10:71]
CB_CT_control_cors <- flattensequareMatrix(cor.prob(CB_CT_dataset[1:6,])) # control cors
CB_CT_transgenic_cors <- flattensequareMatrix(cor.prob(CB_CT_dataset[7:12,])) # transgenic cors
CB_C_and_T_cors <- cbind(CB_CT_control_cors,CB_CT_transgenic_cors[3:4])
colnames(CB_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

CB_CT_fisher = as.array(my_compcorr(CB_C_and_T_cors))
CB_C_and_T_cors_compcor = cbind(CB_C_and_T_cors, "fisher-r-to-z" = CB_CT_fisher)

#############  find all T and C cors for LIV tissue dataset
liv_CT_dataset <- metab_liv[1:12, 10:71]
liv_CT_control_cors <- flattensequareMatrix(cor.prob(liv_CT_dataset[1:6,])) # control cors
liv_CT_transgenic_cors <- flattensequareMatrix(cor.prob(liv_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
liv_C_and_T_cors <- cbind(liv_CT_control_cors,liv_CT_transgenic_cors[3:4])
colnames(liv_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column
liv_CT_fisher = as.array(my_compcorr(liv_C_and_T_cors))
liv_C_and_T_cors_compcor = cbind(liv_C_and_T_cors, "fisher-r-to-z" = liv_CT_fisher)

#############  find all T and C cors for MCTX tissue dataset
mctx_CT_dataset <- metab_mctx[1:12, 10:71]
mctx_CT_control_cors <- flattensequareMatrix(cor.prob(mctx_CT_dataset[1:6,])) # control cors
mctx_CT_transgenic_cors <- flattensequareMatrix(cor.prob(mctx_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
mctx_C_and_T_cors <- cbind(mctx_CT_control_cors,mctx_CT_transgenic_cors[3:4])
colnames(mctx_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column
mctx_CT_fisher = as.array(my_compcorr(mctx_C_and_T_cors))
mctx_C_and_T_cors_compcor = cbind(mctx_C_and_T_cors, "fisher-r-to-z" = mctx_CT_fisher)

#############  find all T and C cors for HIPP tissue dataset
hipp_CT_dataset <- metab_hipp[1:12, 10:70]
hipp_CT_control_cors <- flattensequareMatrix(cor.prob(hipp_CT_dataset[1:6,])) # control cors
hipp_CT_transgenic_cors <- flattensequareMatrix(cor.prob(hipp_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
hipp_C_and_T_cors <- cbind(hipp_CT_control_cors,hipp_CT_transgenic_cors[3:4])
colnames(hipp_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column
hipp_CT_fisher = as.array(my_compcorr(hipp_C_and_T_cors))
hipp_C_and_T_cors_compcor = cbind(hipp_C_and_T_cors, "fisher-r-to-z" = hipp_CT_fisher)

################### Import Biocrates Metabolite data #####################

biocrates_all <- read.csv("EM_5yr_Biocrates_allMetabs_MG.csv", stringsAsFactors = FALSE)

# Separate Biocrates datasets based on tissue type
biocrates_plasma <- biocrates_all[biocrates_all$tissue == "LC-MS Metabolites Plasma", ]; rownames(biocrates_plasma) <- biocrates_plasma$sample
biocrates_mctx <- biocrates_all[biocrates_all$tissue == "LC-MS Metabolites Motor Cortex", ]; rownames(biocrates_mctx) <- biocrates_mctx$sample
biocrates_liv <- biocrates_all[biocrates_all$tissue == "LC-MS Metabolites Liver", ]; rownames(biocrates_liv) <- biocrates_liv$sample
biocrates_cb <- biocrates_all[biocrates_all$tissue == "LC-MS Metabolites Cerebellum", ]; rownames(biocrates_cb) <- biocrates_cb$sample

#############  find all T and C cors for CB tissue dataset
CB_biocrates_CT_dataset <- biocrates_cb[1:12, 6:173]
CB_biocrates_CT_control_cors <- flattensequareMatrix(cor.prob(CB_biocrates_CT_dataset[1:6,])) # control cors
CB_biocrates_CT_transgenic_cors <- flattensequareMatrix(cor.prob(CB_biocrates_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
CB_biocrates_C_and_T_cors <- cbind(CB_biocrates_CT_control_cors, CB_biocrates_CT_transgenic_cors[3:4])
colnames(CB_biocrates_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

CB_biocrates_CT_fisher = as.array(my_compcorr(CB_biocrates_C_and_T_cors))
CB_biocrates_C_and_T_cors_compcor = cbind(CB_biocrates_C_and_T_cors, "fisher-r-to-z" = CB_biocrates_CT_fisher)

#############  find all T and C cors for mctx tissue dataset
mctx_biocrates_CT_dataset <- biocrates_mctx[1:12, 6:173]
mctx_biocrates_CT_control_cors <- flattensequareMatrix(cor.prob(mctx_biocrates_CT_dataset[1:6,])) # control cors
mctx_biocrates_CT_transgenic_cors <- flattensequareMatrix(cor.prob(mctx_biocrates_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
mctx_biocrates_C_and_T_cors <- cbind(mctx_biocrates_CT_control_cors, mctx_biocrates_CT_transgenic_cors[3:4])
colnames(mctx_biocrates_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

mctx_biocrates_CT_fisher = as.array(my_compcorr(mctx_biocrates_C_and_T_cors))
mctx_biocrates_C_and_T_cors_compcor = cbind(mctx_biocrates_C_and_T_cors, "fisher-r-to-z" = mctx_biocrates_CT_fisher)

#############  find all T and C cors for plasma tissue dataset
plasma_biocrates_CT_dataset <- biocrates_plasma[1:12, 6:173]
plasma_biocrates_CT_control_cors <- flattensequareMatrix(cor.prob(plasma_biocrates_CT_dataset[1:6,])) # control cors
plasma_biocrates_CT_transgenic_cors <- flattensequareMatrix(cor.prob(plasma_biocrates_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
plasma_biocrates_C_and_T_cors <- cbind(plasma_biocrates_CT_control_cors, plasma_biocrates_CT_transgenic_cors[3:4])
colnames(plasma_biocrates_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

plasma_biocrates_CT_fisher = as.array(my_compcorr(plasma_biocrates_C_and_T_cors))
plasma_biocrates_C_and_T_cors_compcor = cbind(plasma_biocrates_C_and_T_cors, "fisher-r-to-z" = plasma_biocrates_CT_fisher)

#############  find all T and C cors for liv tissue dataset
liv_biocrates_CT_dataset <- biocrates_liv[1:12, 6:173]
liv_biocrates_CT_control_cors <- flattensequareMatrix(cor.prob(liv_biocrates_CT_dataset[1:6,])) # control cors
liv_biocrates_CT_transgenic_cors <- flattensequareMatrix(cor.prob(liv_biocrates_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
liv_biocrates_C_and_T_cors <- cbind(liv_biocrates_CT_control_cors, liv_biocrates_CT_transgenic_cors[3:4])
colnames(liv_biocrates_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

liv_biocrates_CT_fisher = as.array(my_compcorr(liv_biocrates_C_and_T_cors))
liv_biocrates_C_and_T_cors_compcor = cbind(liv_biocrates_C_and_T_cors, "fisher-r-to-z" = liv_biocrates_CT_fisher)

#######################################################################################################################
#################################################### Urea Data ########################################################
#######################################################################################################################

urea_all_full = read.csv("EM_USE THIS 2012 5yr cohort all Urea related data combined_MG.csv", stringsAsFactors = FALSE)
urea_all = urea_all_full[c(1:72, 85:139),1:9]

# Separating urea data by tissue type
urea_serum <- urea_all[urea_all$tissue == "Urea Serum", ]; rownames(urea_serum) <- urea_serum$sample
urea_serum$tissue = "Serum"
urea_urine <- urea_all[urea_all$tissue == "Urea Urine", ]; rownames(urea_urine) <- urea_urine$sample
urea_urine$tissue = "Urine"
urea_cb <- urea_all[urea_all$tissue == "Urea Cerebellum", ]; rownames(urea_cb) <- urea_cb$sample
urea_cb$tissue = "Cerebellum"
urea_hipp <- urea_all[urea_all$tissue == "Urea Hippocampus", ]; rownames(urea_hipp) <- urea_hipp$sample
urea_hipp$tissue = "Hippocampus"
urea_mctx <- urea_all[urea_all$tissue == "Urea Motor Cortex", ]; rownames(urea_mctx) <- urea_mctx$sample
urea_mctx$tissue = "Cortex"
urea_stri <- urea_all[urea_all$tissue == "Urea Striatum", ]; rownames(urea_stri) <- urea_stri$sample
urea_stri$tissue = "Striatum"
urea_blad <- urea_all[urea_all$tissue == "Urea Bladder", ]; rownames(urea_blad) <- urea_blad$sample
urea_blad$tissue = "Bladder"
urea_heart <- urea_all[urea_all$tissue == "Urea Heart", ]; rownames(urea_heart) <- urea_heart$sample
urea_heart$tissue = "Heart"
urea_kid <- urea_all[urea_all$tissue == "Urea Kidney", ]; rownames(urea_kid) <- urea_kid$sample
urea_kid$tissue = "Kidney"
urea_liv <- urea_all[urea_all$tissue == "Urea Liver", ]; rownames(urea_liv) <- urea_liv$sample
urea_liv$tissue = "Liver"
urea_testes <- urea_all[urea_all$tissue == "Urea Testes", ]; rownames(urea_testes) <- urea_testes$sample
urea_testes$tissue = "Testes"

# Following the same steps taken for previous datasets
urea_serum_CT_dataset <- urea_serum[1:12,8:9]
urea_serum_CT_control_cors <- flattensequareMatrix(cor.prob(urea_serum_CT_dataset[1:6,])) # control cors
urea_serum_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_serum_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_serum_C_and_T_cors <- cbind(urea_serum_CT_control_cors, urea_serum_CT_transgenic_cors[3:4])
colnames(urea_serum_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_serum_CT_fisher = as.array(my_compcorr(urea_serum_C_and_T_cors))
urea_serum_C_and_T_cors_compcor = cbind(urea_serum_C_and_T_cors, "fisher-r-to-z" = urea_serum_CT_fisher)

urea_urine_CT_dataset <- urea_urine[1:12,8:9]
urea_urine_CT_control_cors <- flattensequareMatrix(cor.prob(urea_urine_CT_dataset[1:6,])) # control cors
urea_urine_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_urine_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_urine_C_and_T_cors <- cbind(urea_urine_CT_control_cors, urea_urine_CT_transgenic_cors[3:4])
colnames(urea_urine_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_urine_CT_fisher = as.array(my_compcorr(urea_urine_C_and_T_cors))
urea_urine_C_and_T_cors_compcor = cbind(urea_urine_C_and_T_cors, "fisher-r-to-z" = urea_urine_CT_fisher)

urea_cb_CT_dataset <- urea_cb[1:12,8:9]
urea_cb_CT_control_cors <- flattensequareMatrix(cor.prob(urea_cb_CT_dataset[1:6,])) # control cors
urea_cb_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_cb_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_cb_C_and_T_cors <- cbind(urea_cb_CT_control_cors, urea_cb_CT_transgenic_cors[3:4])
colnames(urea_cb_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_cb_CT_fisher = as.array(my_compcorr(urea_cb_C_and_T_cors))
urea_cb_C_and_T_cors_compcor = cbind(urea_cb_C_and_T_cors, "fisher-r-to-z" = urea_cb_CT_fisher)

urea_hipp_CT_dataset <- urea_hipp[1:12,8:9]
urea_hipp_CT_control_cors <- flattensequareMatrix(cor.prob(urea_hipp_CT_dataset[1:6,])) # control cors
urea_hipp_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_hipp_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_hipp_C_and_T_cors <- cbind(urea_hipp_CT_control_cors, urea_hipp_CT_transgenic_cors[3:4])
colnames(urea_hipp_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_hipp_CT_fisher = as.array(my_compcorr(urea_hipp_C_and_T_cors))
urea_hipp_C_and_T_cors_compcor = cbind(urea_hipp_C_and_T_cors, "fisher-r-to-z" = urea_hipp_CT_fisher)

urea_mctx_CT_dataset <- urea_mctx[1:12,8:9]
urea_mctx_CT_control_cors <- flattensequareMatrix(cor.prob(urea_mctx_CT_dataset[1:6,])) # control cors
urea_mctx_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_mctx_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_mctx_C_and_T_cors <- cbind(urea_mctx_CT_control_cors, urea_mctx_CT_transgenic_cors[3:4])
colnames(urea_mctx_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_mctx_CT_fisher = as.array(my_compcorr(urea_mctx_C_and_T_cors))
urea_mctx_C_and_T_cors_compcor = cbind(urea_mctx_C_and_T_cors, "fisher-r-to-z" = urea_mctx_CT_fisher)

urea_stri_CT_dataset <- urea_stri[1:12,8:9]
urea_stri_CT_control_cors <- flattensequareMatrix(cor.prob(urea_stri_CT_dataset[1:6,])) # control cors
urea_stri_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_stri_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_stri_C_and_T_cors <- cbind(urea_stri_CT_control_cors, urea_stri_CT_transgenic_cors[3:4])
colnames(urea_stri_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_stri_CT_fisher = as.array(my_compcorr(urea_stri_C_and_T_cors))
urea_stri_C_and_T_cors_compcor = cbind(urea_stri_C_and_T_cors, "fisher-r-to-z" = urea_stri_CT_fisher)

urea_blad_CT_dataset <- urea_blad[1:12,8:9]
urea_blad_CT_control_cors <- flattensequareMatrix(cor.prob(urea_blad_CT_dataset[1:6,])) # control cors
urea_blad_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_blad_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_blad_C_and_T_cors <- cbind(urea_blad_CT_control_cors, urea_blad_CT_transgenic_cors[3:4])
colnames(urea_blad_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_blad_CT_fisher = as.array(my_compcorr(urea_blad_C_and_T_cors))
urea_blad_C_and_T_cors_compcor = cbind(urea_blad_C_and_T_cors, "fisher-r-to-z" = urea_blad_CT_fisher)

urea_heart_CT_dataset <- urea_heart[1:12,8:9]
urea_heart_CT_control_cors <- flattensequareMatrix(cor.prob(urea_heart_CT_dataset[1:6,])) # control cors
urea_heart_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_heart_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_heart_C_and_T_cors <- cbind(urea_heart_CT_control_cors, urea_heart_CT_transgenic_cors[3:4])
colnames(urea_heart_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_heart_CT_fisher = as.array(my_compcorr(urea_heart_C_and_T_cors))
urea_heart_C_and_T_cors_compcor = cbind(urea_heart_C_and_T_cors, "fisher-r-to-z" = urea_heart_CT_fisher)

urea_kid_CT_dataset <- urea_kid[1:12,8:9]
urea_kid_CT_control_cors <- flattensequareMatrix(cor.prob(urea_kid_CT_dataset[1:6,])) # control cors
urea_kid_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_kid_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_kid_C_and_T_cors <- cbind(urea_kid_CT_control_cors, urea_kid_CT_transgenic_cors[3:4])
colnames(urea_kid_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_kid_CT_fisher = as.array(my_compcorr(urea_kid_C_and_T_cors))
urea_kid_C_and_T_cors_compcor = cbind(urea_kid_C_and_T_cors, "fisher-r-to-z" = urea_kid_CT_fisher)

urea_liv_CT_dataset <- urea_liv[1:12,8:9]
urea_liv_CT_control_cors <- flattensequareMatrix(cor.prob(urea_liv_CT_dataset[1:6,])) # control cors
urea_liv_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_liv_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
urea_liv_C_and_T_cors <- cbind(urea_liv_CT_control_cors, urea_liv_CT_transgenic_cors[3:4])
colnames(urea_liv_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_liv_CT_fisher = as.array(my_compcorr(urea_liv_C_and_T_cors))
urea_liv_C_and_T_cors_compcor = cbind(urea_liv_C_and_T_cors, "fisher-r-to-z" = urea_liv_CT_fisher)

urea_testes_CT_dataset <- urea_testes[1:7,8:9]
urea_testes_CT_control_cors <- flattensequareMatrix(cor.prob(urea_testes_CT_dataset[1:4,])) # control cors
urea_testes_CT_transgenic_cors <- flattensequareMatrix(cor.prob(urea_testes_CT_dataset[5:7,])) # transgenic cors
# combine 2 dataframes
urea_testes_C_and_T_cors <- cbind(urea_testes_CT_control_cors, urea_testes_CT_transgenic_cors[3:4])
colnames(urea_testes_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

urea_testes_CT_fisher = as.array(my_compcorr(urea_testes_C_and_T_cors))
urea_testes_C_and_T_cors_compcor = cbind(urea_testes_C_and_T_cors, "fisher-r-to-z" = urea_testes_CT_fisher)


##################################################################################################################################################
################################################################ QPCR data #######################################################################
##################################################################################################################################################

qpcr_all <- read.csv("EM_5yrQPCR_allAvg_MG.csv", stringsAsFactors = FALSE)

# Split QPCR data by tissue type
qpcr_mctx <- qpcr_all[qpcr_all$tissue == "QPCR Motor Cortex", ]; rownames(qpcr_mctx) <- qpcr_mctx$sample
qpcr_antistr <- qpcr_all[qpcr_all$tissue == "QPCR Anti-striatum", ]; rownames(qpcr_antistr) <- qpcr_antistr$sample
qpcr_fp <- qpcr_all[qpcr_all$tissue == "QPCR Floor Plate", ]; rownames(qpcr_fp) <- qpcr_fp$sample

#############  find all T and C cors for Motor Cortex QPCR tissue dataset
qpcr_mctx_CT_dataset <- qpcr_mctx[1:12, 10:19]
qpcr_mctx_CT_control_cors <- flattensequareMatrix(cor.prob(qpcr_mctx_CT_dataset[1:6,])) # control cors
qpcr_mctx_CT_transgenic_cors <- flattensequareMatrix(cor.prob(qpcr_mctx_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
qpcr_mctx_C_and_T_cors <- cbind(qpcr_mctx_CT_control_cors, qpcr_mctx_CT_transgenic_cors[3:4])
colnames(qpcr_mctx_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

qpcr_mctx_CT_fisher = as.array(my_compcorr(qpcr_mctx_C_and_T_cors))
qpcr_mctx_C_and_T_cors_compcor = cbind(qpcr_mctx_C_and_T_cors, "fisher-r-to-z" = qpcr_mctx_CT_fisher)

#############  find all T and C cors for QPCR Anti-Striatum tissue dataset
qpcr_antistr_CT_dataset <- qpcr_antistr[1:12, 10:19]
qpcr_antistr_CT_control_cors <- flattensequareMatrix(cor.prob(qpcr_antistr_CT_dataset[1:6,])) # control cors
qpcr_antistr_CT_transgenic_cors <- flattensequareMatrix(cor.prob(qpcr_antistr_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
qpcr_antistr_C_and_T_cors <- cbind(qpcr_antistr_CT_control_cors, qpcr_antistr_CT_transgenic_cors[3:4])
colnames(qpcr_antistr_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

qpcr_antistr_CT_fisher = as.array(my_compcorr(qpcr_antistr_C_and_T_cors))
qpcr_antistr_C_and_T_cors_compcor = cbind(qpcr_antistr_C_and_T_cors, "fisher-r-to-z" = qpcr_antistr_CT_fisher)

#############  find all T and C cors for QPCR Foot Plate tissue dataset
qpcr_fp_CT_dataset <- qpcr_fp[1:12, 10:19]
qpcr_fp_CT_control_cors <- flattensequareMatrix(cor.prob(qpcr_fp_CT_dataset[1:6,])) # control cors
qpcr_fp_CT_transgenic_cors <- flattensequareMatrix(cor.prob(qpcr_fp_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
qpcr_fp_C_and_T_cors <- cbind(qpcr_fp_CT_control_cors, qpcr_fp_CT_transgenic_cors[3:4])
colnames(qpcr_fp_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

qpcr_fp_CT_fisher = as.array(my_compcorr(qpcr_fp_C_and_T_cors))
qpcr_fp_C_and_T_cors_compcor = cbind(qpcr_fp_C_and_T_cors, "fisher-r-to-z" = qpcr_fp_CT_fisher)

##################################################################################################################################################
####################################################### Import proteomics data ###################################################################
##################################################################################################################################################
# These datasets have been separated out by tissue type in excel and do not need to be split in R
# Cerebellum proteomics data
proteo_cb_sig_only = read.csv("Stefano Proteomics Cerebellum Cleaned Significant only Ttest C vs T.csv", stringsAsFactors = FALSE)
rownames(proteo_cb_sig_only) = proteo_cb_sig_only$X
proteo_cb_sig_only$tissue = "Cerebellum Proteomics"
proteo_cb_sig_only_CT_dataset <- proteo_cb_sig_only[1:12, 5:25]
proteo_cb_sig_only_CT_control_cors <- flattensequareMatrix(cor.prob(proteo_cb_sig_only_CT_dataset[1:6,])) # control cors
proteo_cb_sig_only_CT_transgenic_cors <- flattensequareMatrix(cor.prob(proteo_cb_sig_only_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
proteo_cb_sig_only_C_and_T_cors <- cbind(proteo_cb_sig_only_CT_control_cors, proteo_cb_sig_only_CT_transgenic_cors[3:4])
colnames(proteo_cb_sig_only_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

proteo_cb_sig_only_CT_fisher = as.array(my_compcorr(proteo_cb_sig_only_C_and_T_cors))
proteo_cb_sig_only_C_and_T_cors_compcor = cbind(proteo_cb_sig_only_C_and_T_cors, "fisher-r-to-z" = proteo_cb_sig_only_CT_fisher)

##################################################################################################################################################################
# Motor Cortex proteomics dataset
proteo_mctx_sig_only = read.csv("Stefano Proteomics Motor Cortex Cleaned Significant only Ttest C vs T.csv", stringsAsFactors = FALSE)
rownames(proteo_mctx_sig_only) = proteo_mctx_sig_only$X
proteo_mctx_sig_only_CT_dataset <- proteo_mctx_sig_only[1:12, 5:21]
proteo_mctx_sig_only_CT_control_cors <- flattensequareMatrix(cor.prob(proteo_mctx_sig_only_CT_dataset[1:6,])) # control cors
proteo_mctx_sig_only_CT_transgenic_cors <- flattensequareMatrix(cor.prob(proteo_mctx_sig_only_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
proteo_mctx_sig_only_C_and_T_cors <- cbind(proteo_mctx_sig_only_CT_control_cors, proteo_mctx_sig_only_CT_transgenic_cors[3:4])
colnames(proteo_mctx_sig_only_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

proteo_mctx_sig_only_CT_fisher = as.array(my_compcorr(proteo_mctx_sig_only_C_and_T_cors))
proteo_mctx_sig_only_C_and_T_cors_compcor = cbind(proteo_mctx_sig_only_C_and_T_cors, "fisher-r-to-z" = proteo_mctx_sig_only_CT_fisher)

##################################################################################################################################################################
# Striatum proteomics dataset
proteo_str_sig_only = read.csv("Stefano Proteomics Striatum Cleaned Significant only Ttest C vs T.csv", stringsAsFactors = FALSE)
rownames(proteo_str_sig_only) = proteo_str_sig_only$X
proteo_str_sig_only_CT_dataset <- proteo_str_sig_only[1:12, 5:43]
proteo_str_sig_only_CT_control_cors <- flattensequareMatrix(cor.prob(proteo_str_sig_only_CT_dataset[1:6,])) # control cors
proteo_str_sig_only_CT_transgenic_cors <- flattensequareMatrix(cor.prob(proteo_str_sig_only_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
proteo_str_sig_only_C_and_T_cors <- cbind(proteo_str_sig_only_CT_control_cors, proteo_str_sig_only_CT_transgenic_cors[3:4])
colnames(proteo_str_sig_only_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

proteo_str_sig_only_CT_fisher = as.array(my_compcorr(proteo_str_sig_only_C_and_T_cors))
proteo_str_sig_only_C_and_T_cors_compcor = cbind(proteo_str_sig_only_C_and_T_cors, "fisher-r-to-z" = proteo_str_sig_only_CT_fisher)

##################################################################################################################################################################
############################################################# Paul RNASeq data (P-value significant) #############################################################
##################################################################################################################################################################

paul_rnaseq_pval = read.csv("paul_rnaseq_pval_MG.csv", stringsAsFactors = FALSE)
# Setting the rownames of the paul_rnaseq_pval dataset to be the same as the nano_24_DL dataset; they are in the same orientation
rownames(paul_rnaseq_pval) = rownames(nano_24_DL)

paul_rnaseq_pval_CT_dataset <- paul_rnaseq_pval[1:12, 3:820]
paul_rnaseq_pval_CT_control_cors <- flattensequareMatrix(cor.prob(paul_rnaseq_pval_CT_dataset[1:6,])) # control cors
paul_rnaseq_pval_CT_transgenic_cors <- flattensequareMatrix(cor.prob(paul_rnaseq_pval_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
paul_rnaseq_pval_C_and_T_cors <- cbind(paul_rnaseq_pval_CT_control_cors, paul_rnaseq_pval_CT_transgenic_cors[3:4])
colnames(paul_rnaseq_pval_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

paul_rnaseq_pval_CT_fisher = as.array(my_compcorr(paul_rnaseq_pval_C_and_T_cors))
paul_rnaseq_pval_C_and_T_cors_compcor = cbind(paul_rnaseq_pval_C_and_T_cors, "fisher-r-to-z" = paul_rnaseq_pval_CT_fisher)

##################################################################################################################################################################
################################################################ Neuron RNASeq data ##############################################################################
##################################################################################################################################################################

neuron_rnaseq = read.csv("neuron_rnaseq_sig_MG.csv", stringsAsFactors = FALSE)
# Setting the rownames of the neuron_rnaseq dataset to be the same as the nano_24_DL dataset; they are in the same orientation
rownames(neuron_rnaseq) = rownames(nano_24_DL)

neuron_rnaseq_CT_dataset <- neuron_rnaseq[1:12, 3:47]
neuron_rnaseq_CT_control_cors <- flattensequareMatrix(cor.prob(neuron_rnaseq_CT_dataset[1:6,])) # control cors
neuron_rnaseq_CT_transgenic_cors <- flattensequareMatrix(cor.prob(neuron_rnaseq_CT_dataset[7:12,])) # transgenic cors
# combine 2 dataframes
neuron_rnaseq_C_and_T_cors <- cbind(neuron_rnaseq_CT_control_cors, neuron_rnaseq_CT_transgenic_cors[3:4])
colnames(neuron_rnaseq_C_and_T_cors) <- c("Var1", "Var2", "control_cor", "control_pvalue", "transgenic_cor", "transgenic_pvalue")

# then append Fisher value to dataframe as a 7th column

neuron_rnaseq_CT_fisher = as.array(my_compcorr(neuron_rnaseq_C_and_T_cors))
neuron_rnaseq_C_and_T_cors_compcor = cbind(neuron_rnaseq_C_and_T_cors, "fisher-r-to-z" = neuron_rnaseq_CT_fisher)


############## Read in all bootstrap data frames #################
# These files were produced from the R-script 'Bootstrap and permutation app plot code.R'
# The bootstrap data frames for the control and transgenic samples are combined into a single data frame for clean access via the select input variables of the app
nano_24_DM_control_bs.df = read.csv("nano_24_DM_control_bs_df.csv")
nano_24_DM_transgenic_bs.df = read.csv("nano_24_DM_transgenic_bs_df.csv")
nano_24_DM_bs.df = data.frame(cbind(nano_24_DM_control_bs.df, nano_24_DM_transgenic_bs.df))

nano_24_DL_control_bs.df = read.csv("nano_24_DL_control_bs_df.csv")
nano_24_DL_transgenic_bs.df = read.csv("nano_24_DL_transgenic_bs_df.csv")
nano_24_DL_bs.df = data.frame(cbind(nano_24_DL_control_bs.df, nano_24_DL_transgenic_bs.df))

metab_cb_transgenic_bs.df = read.csv("metab_cb_transgenic_bs_df.csv")
metab_cb_control_bs.df = read.csv("metab_cb_control_bs_df.csv")
metab_cb_bs.df = data.frame(cbind(metab_cb_control_bs.df, metab_cb_transgenic_bs.df))

metab_mctx_transgenic_bs.df = read.csv("metab_mctx_transgenic_bs_df.csv")
metab_mctx_control_bs.df = read.csv("metab_mctx_control_bs_df.csv")
metab_mctx_bs.df = data.frame(cbind(metab_mctx_control_bs.df, metab_mctx_transgenic_bs.df))

metab_hipp_transgenic_bs.df = read.csv("metab_hipp_transgenic_bs_df.csv")
metab_hipp_control_bs.df = read.csv("metab_hipp_control_bs_df.csv")
metab_hipp_bs.df = data.frame(cbind(metab_hipp_control_bs.df, metab_hipp_transgenic_bs.df))

metab_liv_transgenic_bs.df = read.csv("metab_liv_transgenic_bs_df.csv")
metab_liv_control_bs.df = read.csv("metab_liv_control_bs_df.csv")
metab_liv_bs.df = data.frame(cbind(metab_liv_control_bs.df, metab_liv_transgenic_bs.df))

biocrates_plasma_transgenic_bs.df = read.csv("biocrates_plasma_transgenic_bs_df.csv")
biocrates_plasma_control_bs.df = read.csv("biocrates_plasma_control_bs_df.csv")
biocrates_plasma_bs.df = data.frame(cbind(biocrates_plasma_control_bs.df, biocrates_plasma_transgenic_bs.df))

biocrates_mctx_transgenic_bs.df = read.csv("biocrates_mctx_transgenic_bs_df.csv")
biocrates_mctx_control_bs.df = read.csv("biocrates_mctx_control_bs_df.csv")
biocrates_mctx_bs.df = data.frame(cbind(biocrates_mctx_control_bs.df, biocrates_mctx_transgenic_bs.df))

biocrates_liv_transgenic_bs.df = read.csv("biocrates_liv_transgenic_bs_df.csv")
biocrates_liv_control_bs.df = read.csv("biocrates_liv_control_bs_df.csv")
biocrates_liv_bs.df = data.frame(cbind(biocrates_liv_control_bs.df, biocrates_liv_transgenic_bs.df))

biocrates_cb_transgenic_bs.df = read.csv("biocrates_cb_transgenic_bs_df.csv")
biocrates_cb_control_bs.df = read.csv("biocrates_cb_control_bs_df.csv")
biocrates_cb_bs.df = data.frame(cbind(biocrates_cb_control_bs.df, biocrates_cb_transgenic_bs.df))

urea_serum_transgenic_bs.df = read.csv("urea_serum_transgenic_bs_df.csv")
urea_serum_control_bs.df = read.csv("urea_serum_control_bs_df.csv")
urea_serum_bs.df = data.frame(cbind(urea_serum_control_bs.df, urea_serum_transgenic_bs.df))

urea_urine_transgenic_bs.df = read.csv("urea_urine_transgenic_bs_df.csv")
urea_urine_control_bs.df = read.csv("urea_urine_control_bs_df.csv")
urea_urine_bs.df = data.frame(cbind(urea_urine_control_bs.df, urea_urine_transgenic_bs.df))

urea_cb_transgenic_bs.df = read.csv("urea_cb_transgenic_bs_df.csv")
urea_cb_control_bs.df = read.csv("urea_cb_control_bs_df.csv")
urea_cb_bs.df = data.frame(cbind(urea_cb_control_bs.df, urea_cb_transgenic_bs.df))

urea_hipp_transgenic_bs.df = read.csv("urea_hipp_transgenic_bs_df.csv")
urea_hipp_control_bs.df = read.csv("urea_hipp_control_bs_df.csv")
urea_hipp_bs.df = data.frame(cbind(urea_hipp_control_bs.df, urea_hipp_transgenic_bs.df))

urea_mctx_transgenic_bs.df = read.csv("urea_mctx_transgenic_bs_df.csv")
urea_mctx_control_bs.df = read.csv("urea_mctx_control_bs_df.csv")
urea_mctx_bs.df = data.frame(cbind(urea_mctx_control_bs.df, urea_mctx_transgenic_bs.df))

urea_stri_transgenic_bs.df = read.csv("urea_stri_transgenic_bs_df.csv")
urea_stri_control_bs.df = read.csv("urea_stri_control_bs_df.csv")
urea_stri_bs.df = data.frame(cbind(urea_stri_control_bs.df, urea_stri_transgenic_bs.df))

urea_blad_transgenic_bs.df = read.csv("urea_blad_transgenic_bs_df.csv")
urea_blad_control_bs.df = read.csv("urea_blad_control_bs_df.csv")
urea_blad_bs.df = data.frame(cbind(urea_blad_control_bs.df, urea_blad_transgenic_bs.df))

urea_heart_transgenic_bs.df = read.csv("urea_heart_transgenic_bs_df.csv")
urea_heart_control_bs.df = read.csv("urea_heart_control_bs_df.csv")
urea_heart_bs.df = data.frame(cbind(urea_heart_control_bs.df, urea_heart_transgenic_bs.df))

urea_kid_transgenic_bs.df = read.csv("urea_kid_transgenic_bs_df.csv")
urea_kid_control_bs.df = read.csv("urea_kid_control_bs_df.csv")
urea_kid_bs.df = data.frame(cbind(urea_kid_control_bs.df, urea_kid_transgenic_bs.df))

urea_liv_transgenic_bs.df = read.csv("urea_liv_transgenic_bs_df.csv")
urea_liv_control_bs.df = read.csv("urea_liv_control_bs_df.csv")
urea_liv_bs.df = data.frame(cbind(urea_liv_control_bs.df, urea_liv_transgenic_bs.df))

urea_testes_transgenic_bs.df = read.csv("urea_testes_transgenic_bs_df.csv")
urea_testes_control_bs.df = read.csv("urea_testes_control_bs_df.csv")
urea_testes_bs.df = data.frame(cbind(urea_testes_control_bs.df, urea_testes_transgenic_bs.df))

qpcr_mctx_transgenic_bs.df = read.csv("qpcr_mctx_transgenic_bs_df.csv")
qpcr_mctx_control_bs.df = read.csv("qpcr_mctx_control_bs_df.csv")
qpcr_mctx_bs.df = data.frame(cbind(qpcr_mctx_control_bs.df, qpcr_mctx_transgenic_bs.df))

qpcr_antistr_transgenic_bs.df = read.csv("qpcr_antistr_transgenic_bs_df.csv")
qpcr_antistr_control_bs.df = read.csv("qpcr_antistr_control_bs_df.csv")
qpcr_antistr_bs.df = data.frame(cbind(qpcr_antistr_control_bs.df, qpcr_antistr_transgenic_bs.df))

qpcr_fp_transgenic_bs.df = read.csv("qpcr_fp_transgenic_bs_df.csv")
qpcr_fp_control_bs.df = read.csv("qpcr_fp_control_bs_df.csv")
qpcr_fp_bs.df = data.frame(cbind(qpcr_fp_control_bs.df, qpcr_fp_transgenic_bs.df))

proteo_cb_sig_only_transgenic_bs.df = read.csv("proteo_cb_sig_only_transgenic_bs_df.csv")
proteo_cb_sig_only_control_bs.df = read.csv("proteo_cb_sig_only_control_bs_df.csv")
proteo_cb_sig_only_bs.df = data.frame(cbind(proteo_cb_sig_only_control_bs.df, proteo_cb_sig_only_transgenic_bs.df))

proteo_mctx_sig_only_transgenic_bs.df = read.csv("proteo_mctx_sig_only_transgenic_bs_df.csv")
proteo_mctx_sig_only_control_bs.df = read.csv("proteo_mctx_sig_only_control_bs_df.csv")
proteo_mctx_sig_only_bs.df = data.frame(cbind(proteo_mctx_sig_only_control_bs.df, proteo_mctx_sig_only_transgenic_bs.df))

proteo_str_sig_only_transgenic_bs.df = read.csv("proteo_str_sig_only_transgenic_bs_df.csv")
proteo_str_sig_only_control_bs.df = read.csv("proteo_str_sig_only_control_bs_df.csv")
proteo_str_sig_only_bs.df = data.frame(cbind(proteo_str_sig_only_control_bs.df, proteo_str_sig_only_transgenic_bs.df))

paul_rnaseq_pval_transgenic_bs.df = read.csv("paul_rnaseq_pval_transgenic_bs_df.csv")
paul_rnaseq_pval_control_bs.df = read.csv("paul_rnaseq_pval_control_bs_df.csv")
paul_rnaseq_pval_bs.df = data.frame(cbind(paul_rnaseq_pval_control_bs.df, paul_rnaseq_pval_transgenic_bs.df))

#liver_rnaseq_transgenic_bs.df = read.csv("liver_rnaseq_transgenic_bs_df.csv")
#liver_rnaseq_control_bs.df = read.csv("liver_rnaseq_control_bs_df.csv")
#liver_rnaseq_bs.df = data.frame(cbind(liver_rnaseq_control_bs.df, liver_rnaseq_transgenic_bs.df))

neuron_rnaseq_transgenic_bs.df = read.csv("neuron_rnaseq_transgenic_bs_df.csv")
neuron_rnaseq_control_bs.df = read.csv("neuron_rnaseq_control_bs_df.csv")
neuron_rnaseq_bs.df = data.frame(cbind(neuron_rnaseq_control_bs.df, neuron_rnaseq_transgenic_bs.df))


######## Read in Permutation data frames #########

nano_24_DL_perm.df = read.csv("nano_24_DL_perm_df.csv")
nano_24_DM_perm.df = read.csv("nano_24_DM_perm_df.csv")
metab_cb_perm.df = read.csv("metab_cb_perm_df.csv")
metab_mctx_perm.df = read.csv("metab_mctx_perm_df.csv")
metab_liv_perm.df = read.csv("metab_liv_perm_df.csv")
metab_hipp_perm.df = read.csv("metab_hipp_perm_df.csv")
biocrates_cb_perm.df = read.csv("biocrates_cb_perm_df.csv")
biocrates_liv_perm.df = read.csv("biocrates_liv_perm_df.csv")
biocrates_mctx_perm.df = read.csv("biocrates_mctx_perm_df.csv")
biocrates_plasma_perm.df = read.csv("biocrates_plasma_perm_df.csv")
#qpcr_mctx_perm.df = read.csv("qpcr_mctx_perm_df.csv")
#qpcr_antistr_perm.df = read.csv("qpcr_antistr_perm_df.csv")
#qpcr_fp_perm.df = read.csv("qpcr_fp_perm_df.csv")
urea_serum_perm.df = read.csv("urea_serum_perm_df.csv")
urea_urine_perm.df = read.csv("urea_urine_perm_df.csv")
urea_cb_perm.df = read.csv("urea_cb_perm_df.csv")
urea_hipp_perm.df = read.csv("urea_hipp_perm_df.csv")
urea_mctx_perm.df = read.csv("urea_mctx_perm_df.csv")
urea_stri_perm.df = read.csv("urea_stri_perm_df.csv")
urea_blad_perm.df = read.csv("urea_blad_perm_df.csv")
urea_heart_perm.df = read.csv("urea_heart_perm_df.csv")
urea_kid_perm.df = read.csv("urea_kid_perm_df.csv")
urea_liv_perm.df = read.csv("urea_liv_perm_df.csv")
urea_testes_perm.df = read.csv("urea_testes_perm_df.csv")
proteo_cb_sig_only_perm.df = read.csv("proteo_cb_sig_only_perm_df.csv")
proteo_mctx_sig_only_perm.df = read.csv("proteo_mctx_sig_only_perm_df.csv")
proteo_str_sig_only_perm.df = read.csv("proteo_str_sig_only_perm_df.csv")
paul_rnaseq_pval_perm.df = read.csv("paul_rnaseq_pval_perm_df.csv")
#liver_rnaseq_perm.df = read.csv("liver_rnaseq_perm_df.csv")
neuron_rnaseq_perm.df = read.csv("neuron_rnaseq_perm_df.csv")

##################################################################################################################################################################
# This list is used when generating the PCA plot and the arrows to represent the selected variables
# The variables that do not have data collected for them are removed from the datasets to prevent errors in the PCA calculations
data_sets_list = list("Dorsolateral" = nano_24_DL[,5:30], 
                      "Dorsomedial" = nano_24_DM[,5:30],
                      "Metabolites (GC-MS) Cerebellum" = metab_cb[,10:71][,-which(colSums(metab_cb[,10:71]) == 0)],
                      "Metabolites (GC-MS) Motor Cortex" = metab_mctx[,10:71][,-which(colSums(metab_mctx[,10:71]) == 0)],
                      "Metabolites (GC-MS) Liver" = metab_liv[,10:71][,-which(colSums(metab_liv[,10:71]) == 0)],
                      "Metabolites (GC-MS) Hippocampus" = metab_hipp[,10:70][,-which(colSums(metab_hipp[,10:70]) == 0)],
                      "Metabolites (LC-MS) Cerebellum" = biocrates_cb[,6:173][,-which(colSums(biocrates_cb[,6:173], na.rm = T) == 0 | is.na(colSums(biocrates_cb[,6:173])))],
                      "Metabolites (LC-MS) Motor Cortex" = biocrates_mctx[,6:173][,-which(colSums(biocrates_mctx[,6:173], na.rm = T) == 0 | is.na(colSums(biocrates_mctx[,6:173])))],
                      "Metabolites (LC-MS) Liver" = biocrates_liv[,6:173][,-which(colSums(biocrates_liv[,6:173], na.rm = T) == 0 | is.na(colSums(biocrates_liv[,6:173])))],
                      "Metabolites (LC-MS) Plasma" = biocrates_plasma[,6:173][,-which(colSums(biocrates_plasma[,6:173], na.rm = T) == 0 | is.na(colSums(biocrates_plasma[,6:173])))],
                      #"qPCR Motor Cortex" = qpcr_mctx[,10:19][,-which(colSums(qpcr_mctx[,10:19]) == 0)],
                      #"qPCR Anti-striatum" = qpcr_antistr[,10:19][,-which(colSums(qpcr_antistr[,10:19]) == 0)],
                      #"qPCR Foot Plate" = qpcr_fp[,10:19][,-which(colSums(qpcr_fp[,10:19]) == 0)],
                      "Urea Serum" = urea_serum[,8:9],
                      "Urea Urine" = urea_urine[,8:9],
                      "Urea Cerebellum" = urea_cb[,8:9][,-which(colSums(urea_cb[,8:9]) == 0 | is.na(colSums(urea_cb[,8:9])))],
                      "Urea Hippocampus" = urea_hipp[,8:9],
                      "Urea Motor Cortex" = urea_mctx[,8:9][,-which(colSums(urea_mctx[,8:9]) == 0 | is.na(colSums(urea_mctx[,8:9])))],
                      "Urea Striatum" = urea_stri[,8:9],
                      "Urea Bladder" = urea_blad[,8:9],
                      "Urea Heart" = urea_heart[,8:9][,-which(colSums(urea_heart[,8:9]) == 0 | is.na(colSums(urea_heart[,8:9])))],
                      "Urea Kidney" = urea_kid[,8:9],
                      "Urea Liver" = urea_liv[,8:9][,-which(colSums(urea_liv[,8:9]) == 0 | is.na(colSums(urea_liv[,8:9])))],
                      "Urea Testes" = urea_testes[,8:9],
                      "Cerebellum Proteomics" = proteo_cb_sig_only[,4:24],
                      "Motor Cortex Proteomics" = proteo_mctx_sig_only[,4:20],
                      "Striatum Proteomics" = proteo_str_sig_only[,4:42],
                      "Striatum RNASeq data" = paul_rnaseq_pval[,3:819],
                      #"Liver RNASeq" = liver_rnaseq[,3:67],
                      "Striatal Neuron RNASeq data" = neuron_rnaseq[,3:47]
)

# This list is used when producing the bootstrap plots in order to produce the vline and mean value
data_sets_list2 = list("nano_24_DL_bs.df" = nano_24_DL[,5:30], 
                       "nano_24_DM_bs.df" = nano_24_DM[,5:30],
                       "metab_cb_bs.df" = metab_cb[,10:71],
                       "metab_mctx_bs.df" = metab_mctx[,10:71],
                       "metab_liv_bs.df" = metab_liv[,10:71],
                       "metab_hipp_bs.df" = metab_hipp[,10:70],
                       "biocrates_cb_bs.df" = biocrates_cb[,6:173],
                       "biocrates_mctx_bs.df" = biocrates_mctx[,6:173],
                       "biocrates_liv_bs.df" = biocrates_liv[,6:173],
                       "biocrates_plasma_bs.df" = biocrates_plasma[,6:173],
                       #"qpcr_mctx_bs.df" = qpcr_mctx[,10:19],
                       #"qpcr_antistr_bs.df" = qpcr_antistr[,10:19],
                       #"qpcr_fp_bs.df" = qpcr_fp[,10:19],
                       "urea_serum_bs.df" = urea_serum[,8:9],
                       "urea_urine_bs.df" = urea_urine[,8:9],
                       "urea_cb_bs.df" = urea_cb[,8:9],
                       "urea_hipp_bs.df" = urea_hipp[,8:9],
                       "urea_mctx_bs.df" = urea_mctx[,8:9],
                       "urea_stri_bs.df" = urea_stri[,8:9],
                       "urea_blad_bs.df" = urea_blad[,8:9],
                       "urea_heart_bs.df" = urea_heart[,8:9],
                       "urea_kid_bs.df" = urea_kid[,8:9],
                       "urea_liv_bs.df" = urea_liv[,8:9],
                       "urea_testes_bs.df" = urea_testes[,8:9],
                       "proteo_cb_sig_only_bs.df" = proteo_cb_sig_only[,4:24],
                       "proteo_mctx_sig_only_bs.df" = proteo_mctx_sig_only[,4:20],
                       "proteo_str_sig_only_bs.df" = proteo_str_sig_only[,4:42],
                       "paul_rnaseq_pval_bs.df" = paul_rnaseq_pval[,3:819],
                       #"liver_rnaseq_bs.df" = liver_rnaseq[,3:67],
                       "neuron_rnaseq_bs.df" = neuron_rnaseq[,3:47]
)

# This list is used when producing the permutation plots in order to produce the vline and significance p-value
data_sets_list3 = list("nano_24_DL_bs.df" = nano_24_DL_perm.df,
                       "nano_24_DM_bs.df" = nano_24_DM_perm.df,
                       "metab_cb_bs.df" = metab_cb_perm.df,
                       "metab_mctx_bs.df" = metab_mctx_perm.df,
                       "metab_liv_bs.df" = metab_liv_perm.df,
                       "metab_hipp_bs.df" = metab_hipp_perm.df,
                       "biocrates_cb_bs.df" = biocrates_cb_perm.df,
                       "biocrates_mctx_bs.df" = biocrates_mctx_perm.df,
                       "biocrates_liv_bs.df" = biocrates_liv_perm.df,
                       "biocrates_plasma_bs.df" = biocrates_plasma_perm.df,
                       #"qpcr_mctx_bs.df" = qpcr_mctx_perm.df,
                       #"qpcr_antistr_bs.df" = qpcr_antistr_perm.df,
                       #"qpcr_fp_bs.df" = qpcr_fp_perm.df,
                       "urea_serum_bs.df" = urea_serum_perm.df,
                       "urea_urine_bs.df" = urea_urine_perm.df,
                       "urea_cb_bs.df" = urea_cb_perm.df,
                       "urea_hipp_bs.df" = urea_hipp_perm.df,
                       "urea_mctx_bs.df" = urea_mctx_perm.df,
                       "urea_stri_bs.df" = urea_stri_perm.df,
                       "urea_blad_bs.df" = urea_blad_perm.df,
                       "urea_heart_bs.df" = urea_heart_perm.df,
                       "urea_kid_bs.df" = urea_kid_perm.df,
                       "urea_liv_bs.df" = urea_liv_perm.df,
                       "urea_testes_bs.df" = urea_testes_perm.df,
                       "proteo_cb_sig_only_bs.df" = proteo_cb_sig_only_perm.df,
                       "proteo_mctx_sig_only_bs.df" = proteo_mctx_sig_only_perm.df,
                       "proteo_str_sig_only_bs.df" = proteo_str_sig_only_perm.df,
                       "paul_rnaseq_pval_bs.df" = paul_rnaseq_pval_perm.df,
                       #"liver_rnaseq_bs.df" = liver_rnaseq_perm.df,
                       "neuron_rnaseq_bs.df" = neuron_rnaseq_perm.df
)

C_and_T_compcor_list = list("Dorsolateral" = DL_C_and_T_cors_compcor,
                            "Dorsomedial" = DM_C_and_T_cors_compcor,
                            "Metabolites (GC-MS) Cerebellum" = CB_C_and_T_cors_compcor,
                            "Metabolites (GC-MS) Hippocampus" = hipp_C_and_T_cors_compcor,
                            "Metabolites (GC-MS) Liver" = liv_C_and_T_cors_compcor,
                            "Metabolites (GC-MS) Motor Cortex" = mctx_C_and_T_cors_compcor,
                            "Metabolites (LC-MS) Motor Cortex" = mctx_biocrates_C_and_T_cors_compcor,
                            "Metabolites (LC-MS) Cerebellum" = CB_biocrates_C_and_T_cors_compcor,
                            "Metabolites (LC-MS) Plasma" = plasma_biocrates_C_and_T_cors_compcor,
                            "Metabolites (LC-MS) Liver" = liv_biocrates_C_and_T_cors_compcor,
                            #"QPCR Motor Cortex" = qpcr_mctx_C_and_T_cors_compcor,
                            #"QPCR Anti-striatum" = qpcr_antistr_C_and_T_cors_compcor,
                            #"QPCR Floor Plate" = qpcr_fp_C_and_T_cors_compcor,
                            "Urea Serum" = urea_serum_C_and_T_cors_compcor,
                            "Urea Urine" = urea_urine_C_and_T_cors_compcor,
                            "Urea Cerebellum" = urea_cb_C_and_T_cors_compcor,
                            "Urea Hippocampus" = urea_hipp_C_and_T_cors_compcor,
                            "Urea Motor Cortex" = urea_mctx_C_and_T_cors_compcor,
                            "Urea Striatum" = urea_stri_C_and_T_cors_compcor,
                            "Urea Bladder" = urea_blad_C_and_T_cors_compcor,
                            "Urea Heart" = urea_heart_C_and_T_cors_compcor,
                            "Urea Kidney" = urea_kid_C_and_T_cors_compcor,
                            "Urea Liver" = urea_liv_C_and_T_cors_compcor,
                            "Urea Testes" = urea_testes_C_and_T_cors_compcor,
                            "Cerebellum Proteomics" = proteo_cb_sig_only_C_and_T_cors_compcor,
                            "Motor Cortex Proteomics" = proteo_mctx_sig_only_C_and_T_cors_compcor,
                            "Striatum Proteomics" = proteo_str_sig_only_C_and_T_cors_compcor,
                            "Striatum RNASeq data" = paul_rnaseq_pval_C_and_T_cors_compcor,
                            #"Liver RNASeq" = liver_rnaseq_C_and_T_cors_compcor,
                            "Striatal Neuron RNASeq" = neuron_rnaseq_C_and_T_cors_compcor
)

data_set_names = list("nano_24" = "nanoString Data", 
                      "metab_all" = "Metabolites (GC-MS)",
                      "biocrates_all" = "Metabolites (LC-MS)",
                      "urea_all" = "Urea Data",
                      #"qpcr_all" = "QPCR Data",
                      "CB_proteomics_download.df" = "Cerebellum Proteomics",
                      "MCTX_proteomics_download.df" = "Motor Cortex Proteomics",
                      "STRI_proteomics_download.df" = "Striatum Proteomics",
                      "paul_rnaseq_pval" = "Striatum RNASeq data",
                      "neuron_rnaseq" = "Striatal Neuron data")

# Read in full proteomics datasets
CB_proteomics_download.df = read.csv("EM_ex stefano 25 6 Proteomics CB 2012 OVT73 (cleaned)_MG.csv")
MCTX_proteomics_download.df = read.csv("EM_ex stefano 25 6 Proteomics MCTX 2012 OVT73 (cleaned)_MG.csv")
STRI_proteomics_download.df = read.csv("EM_ex stefano 25 6 Proteomics STR 2012 OVT73 (cleaned)_MG.csv")

# Read in gene enrichment data frames from DAVID analysis conducted on nanoString, Striatal neuron RNASeq and Striatum RNASeq significant genes. (All from the Striatum list of genes - Paul RNASeq data frame were used as these were previously identified as being significant.)

nanoString_DAVID_analysis_pre <- read.csv("nanoString_DAVID_analysis.csv", stringsAsFactors = FALSE)
Striatum_RNASeq_DAVID_analysis_pre <- read.csv("Striatum_RNASeq_DAVID_analysis.csv", stringsAsFactors = FALSE)
Striatal_Neuron_DAVID_analysis_pre <- read.csv("Striatal_Neuron_DAVID_analysis.csv", stringsAsFactors = FALSE)
# Clean up as the original read in files contain the names of the genes used in the first column, which are subsequently removed.
nanoString_DAVID_analysis = nanoString_DAVID_analysis_pre[1:24, 2:15]
Striatum_RNASeq_DAVID_analysis = Striatum_RNASeq_DAVID_analysis_pre[1:218, 2:15]
Striatal_Neuron_DAVID_analysis = Striatal_Neuron_DAVID_analysis_pre[1:6, 2:15]

#for(i in 1:ncol(nanoString_DAVID_analysis)){
#  if(mode(nanoString_DAVID_analysis[,i]) == 'numeric'){
#    nanoString_DAVID_analysis[,i] = formatC(nanoString_DAVID_analysis[,i], format = "e", digits = 2)
#  }
#}

#for(i in 1:ncol(Striatum_RNASeq_DAVID_analysis)){
#  if(mode(Striatum_RNASeq_DAVID_analysis[,i]) == 'numeric'){
#    Striatum_RNASeq_DAVID_analysis[,i] = formatC(Striatum_RNASeq_DAVID_analysis[,i], format = "e", digits = 2)
#  }
#}

#for(i in 1:ncol(Striatal_Neuron_DAVID_analysis)){
#  if(mode(Striatal_Neuron_DAVID_analysis[,i]) == 'numeric'){
#    Striatal_Neuron_DAVID_analysis[,i] = formatC(Striatal_Neuron_DAVID_analysis[,i], format = "e", digits = 2)
#  }
#}






# Define UI for the application. 
ui <- shinyUI(
  navbarPage(title = "", 
             header = singleton(tags$head(HTML(
               "<script>
                   (function(i,s,o,g,r,a,m){
                   i['GoogleAnalyticsObject']=r;i[r]=i[r]||
                   function(){
                   (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();
                   a=s.createElement(o), m=s.getElementsByTagName(o)[0];
                   a.async=1;
                   a.src=g;m.parentNode.insertBefore(a,m)
                   })
                   (window, document, 'script',
                   '//www.google-analytics.com/analytics.js','ga');

                   ga('create', 'UA-134861665-1', 'auto');
                   ga('send', 'pageview');


                   </script>"
             ),
             tags$style(HTML('.navbar-default {background-color: #00467f;}',
                             '.navbar a{font-weight: bold; color: #fff !important;background-color: #00467f;}',
                             '.navbar-default .navbar-brand{color: #00467f !important;}',
                             '.navbar-nav li a:hover, .navbar-nav > .active > a {color: #00467f !important;
                                               background-color: #FFF !important;
                                               background-image: none !important
                                               }')))),
             
             
             #This is the Home page layout.
                         tabPanel(("Home"),
                                  # HTML code to create the heading line in which the logo will be contained
                                  HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:10px;position:absolute;display:inline-block;font-size:50px;'>HDSheep</h1> <div class='UOAlogos' style='overflow: hidden;display:inline-block;margin-left:40%'> <img src='UOA.png' style='width:211.2px;height:96px;'/><img src='CBR.png' style='width:211.2px;height:96px;'/></div><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  # HTML code to create the background image
                                  #HTML("<style> body {background: url('sheep_bg20.png'); background-attachment: fixed; background-position: center top; background-repeat: no-repeat; width: 100%; background-size:cover;}</style>"),
                                  sidebarLayout(
                                    # Create a side panel
                                    sidebarPanel(
                                      # Text in the side panel
                                      helpText(h3("Links to journals:")),
                                      # Titles and links to relevnt journal articles
                                      helpText(h6("1) miRNA treatment of mutant HTT")),
                                      a("Artificial miRNAs reduce human mutant Huntingtin throughout the striatum in a transgenic sheep model of Huntington's disease (Pfister E, et al. 2017).", href="https://www.ncbi.nlm.nih.gov/pubmed/29207890"),
                                      br(),
                                      helpText(h6("2) Brain urea identified in OVT73 brain tissue:")),
                                      a("Brain urea increase is an early Huntington's disease pathogenic event observed in a prodromal transgenic sheep model and HD cases  (Handley R, et al. 2017)", href="https://www.ncbi.nlm.nih.gov/pubmed/29229845"),
                                      br(),
                                      helpText(h6("3) Metabolic profiling of the OVT73 model:")),
                                      a("Metabolic profiling of presymptomatic Huntington’s disease sheep reveals novel biomarkers (Skene D, et al. 2017)", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320451/"),
                                      br(),
                                      helpText(h6("4) Metabolic disruption identified in OVT73:")),
                                      a("Metabolic disruption identified in the Huntington's disease transgenic sheep model (Handley R, et al. 2016)", href="https://www.ncbi.nlm.nih.gov/pubmed/26864449"),
                                      br(),
                                      helpText(h6("5) Molecular characterisation of the OVT73 model:")),
                                      a("Early and progressive circadian abnormalities in Huntington's disease sheep are unmasked by social environment (Morton A, et al. 2014)", href="https://www.ncbi.nlm.nih.gov/pubmed/24488771"),
                                      br(),
                                      helpText(h6("6) Molecular characterisation of the OVT73 model:")),
                                      a("Further molecular characterisation of the OVT73 transgenic sheep model of Huntington's disease identifies cortical aggregates (Reid S, et al. 2013)", href="https://www.ncbi.nlm.nih.gov/pubmed/25062676"),
                                      br(),
                                      helpText(h6("7) Development of the OVT73 transgenic sheep line:")),
                                      a("An ovine transgenic Huntington's disease model (Jacobsen J, et al. 2010)", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2860888/"),
                                      br(),
                                      br(),
                                      helpText(h6("OMIM entries for Huntington's disease and the Huntingtin gene:")),
                                      a("Huntington's disease OMIM entry", href = "https://www.omim.org/entry/143100?search=Huntington%27s%20disease&highlight=huntington%20disease"),
                                      br(),
                                      a("Huntingtin Gene OMIM entry", href = "https://www.omim.org/entry/613004?search=Huntington%27s%20disease&highlight=huntington%20disease"),
                                      br()
                                    ), position = 'right',
                                    # Create the main panel of the page
                                    mainPanel(
                                      # HTML code to insert the sheep image
                                      HTML('<p><img src="Sheep.png" style="width:500px;height:375px;margin-left:15px;margin-bottom:15px;float: right;"> <h2>Welcome to the HD Sheep Model (OVT73) data exploration tool. </h2> Huntington\'s Disease (HD) is an autosomal dominantly inherited genetic disorder characterised by spontaneous movements, cognitive impairment, psychiatric disturbance and progressive neurodegeneration <a href = "https://www.omim.org/entry/143100"> OMIM # 143100</a>.In collaboration with our partners <br></p><img src="BRNZ.png" style="width:90px;height:90px;display:inline-block;margin-top:25px;margin-left:10px"> <img src="CHDI.png" style="width:162px;height:90px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="SARDI.png" style="width:54px;height:90px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="UniofAdelaide.png" style="width:144px;height:90px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="SAHMRI.png" style="width:162px;height:90px;display:inline-block;margin-top:25px;margin-left:1%">'),
                                      br(),
                                      # Two sets of text input, had to be separated due to layout issues
                                      br()
                                      
                                    ))
                         ),
             
             tabPanel(("Datasets"),
                      HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Datasets Presented in this Project</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                      sidebarLayout(
                        sidebarPanel(
                          helpText(h3("Links to journals:")),
                          helpText(h6("1) miRNA treatment of mutant HTT")),
                          a("Artificial miRNAs reduce human mutant Huntingtin throughout the striatum in a transgenic sheep model of Huntington's disease (Pfister E, et al. 2017).", href="https://www.ncbi.nlm.nih.gov/pubmed/29207890"),
                          br(),
                          helpText(h6("2) Brain urea identified in OVT73 brain tissue:")),
                          a("Brain urea increase is an early Huntington's disease pathogenic event observed in a prodromal transgenic sheep model and HD cases  (Handley R, et al. 2017)", href="https://www.ncbi.nlm.nih.gov/pubmed/29229845"),
                          br(),
                          helpText(h6("3) Metabolic profiling of the OVT73 model:")),
                          a("Metabolic profiling of presymptomatic Huntington’s disease sheep reveals novel biomarkers (Skene D, et al. 2017)", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320451/"),
                          br(),
                          helpText(h6("4) Metabolic disruption identified in OVT73:")),
                          a("Metabolic disruption identified in the Huntington's disease transgenic sheep model (Handley R, et al. 2016)", href="https://www.ncbi.nlm.nih.gov/pubmed/26864449"),
                          br(),
                          helpText(h6("5) Molecular characterisation of the OVT73 model:")),
                          a("Early and progressive circadian abnormalities in Huntington's disease sheep are unmasked by social environment (Morton A, et al. 2014)", href="https://www.ncbi.nlm.nih.gov/pubmed/24488771"),
                          br(),
                          helpText(h6("6) Molecular characterisation of the OVT73 model:")),
                          a("Further molecular characterisation of the OVT73 transgenic sheep model of Huntington's disease identifies cortical aggregates (Reid S, et al. 2013)", href="https://www.ncbi.nlm.nih.gov/pubmed/25062676"),
                          br(),
                          helpText(h6("7) Development of the OVT73 transgenic sheep line:")),
                          a("An ovine transgenic Huntington's disease model (Jacobsen J, et al. 2010)", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2860888/"),
                          br(),
                          br(),
                          helpText(h6("OMIM entries for Huntington's disease and the Huntingtin gene:")),
                          a("Huntington's disease OMIM entry", href = "https://www.omim.org/entry/143100?search=Huntington%27s%20disease&highlight=huntington%20disease"),
                          br(),
                          a("Huntingtin Gene OMIM entry", href = "https://www.omim.org/entry/613004?search=Huntington%27s%20disease&highlight=huntington%20disease"),
                          br()), position = 'right', 
                        mainPanel(HTML('<div style="display:block;width:100%;"><div style="position:relative;right:0px;font-size:14px;">The datasets presented in this application have been collected from a single cohort of 5-year-old OVT73 (n=6) and control (n=6) sheep. Harvest information for each of the 12 samples is provided in the table below. An overview of each dataset, including experimental design, is in preparation.<br><br>Please look at the presented journal papers here and on the home page to read on the research that has previously been conducted on these datasets.<br><br>Individual datasets can be downloaded using the drop-down box below:<br><br></div></div>'),
                                  HTML('<p><img src="Sheep_info.png" style="width:100%;margin-left:0px;"><br><br>'),
                                  helpText(h5("Datasets can be downloaded using the drop-down box below:")),
                                  selectInput("dataset", "",
                                              choices = list("Transcript expression (nanoString)" = "nano_24", 
                                                             #"Transcript expression (QPCR)" = "qpcr_all",
                                                             "Transcript expression (RNAseq Striatum)" = "paul_rnaseq_pval",
                                                             "Transcript expression (RNAseq Striatal Neurons)" = "neuron_rnaseq",
                                                             "Proteomics (Cerebellum)" = "CB_proteomics_download.df",
                                                             "Proteomics (Motor Cortex)" = "MCTX_proteomics_download.df",
                                                             "Proteomics (Striatum)" = "STRI_proteomics_download.df",
                                                             "Metabolites (GC-MS)" = "metab_all",
                                                             "Metabolites (LC-MS)" = "biocrates_all",
                                                             "Urea Quantification" = "urea_all"
                                                             )),
                                  downloadButton("downloadData", "Download Dataset"),
                                  br(),
                                  br(),
                                  br(),
                                  br(),
                                  br(),
                                  br(),
                                  br(),
                                  br(),
                                  br())
                      )),
             
                         # Second tab
                         tabPanel(("Student's T-test"),
                                  HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Student's T-test</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  #Creating the layout of the sidebar on the Single-variable statistics page.
                                  sidebarLayout(
                                    sidebarPanel(#Creates a selection box.
                                                 selectInput("select_1", label = h5("Select dataset"),
                                                             choices = list("Transcript expression (nanoString)" = "nano_24", 
                                                                            #"Transcript expression (QPCR)" = "qpcr_all",
                                                                            "Transcript expression (RNAseq Striatum)" = "paul_rnaseq_pval",
                                                                            "Transcript expression (RNAseq Striatal Neurons)" = "neuron_rnaseq",
                                                                            "Proteomics" = "Proteomics",
                                                                            "Metabolites (GC-MS)" = "metab_all",
                                                                            "Metabolites (LC-MS)" = "biocrates_all",
                                                                            "Urea Quantification" = "urea_all")
                                                 ),
                                                 #Creates a selection box.
                                                 selectInput("select_2", label = h5("Select tissue"), 
                                                             choices = list("Dorsolateral" = "nano_24_DL", 
                                                                            "Dorsomedial" = "nano_24_DM")
                                                 ),
                                                 #Creates a selection box. This takes the names from the dataset columns as the selection variables and the data
                                                 #associated with the variable.
                                                 selectInput("select_3", label = h5("Select specific variable of interest (transcript / metabolite / protein)"), 
                                                             choices = names(nano_24)[6:31]),
                                                 downloadButton("pdflink1")
                                    ),
                                    # This is the layout of the single variable statistics page.
                                    mainPanel(
                                      HTML('<div class="container"; style="display:block;width:100%;"><button {class="more/less1" data-toggle="collapse" data-target="#description1" style="position: relative; left:0px; border-radius: 5px; background-color: #00467f;border: none; color: white; padding: 5px 8px; text-align: center; text-decoration: none;font-size: 12px; margin-bottom: 10px; cursor: pointer;"  >More/Less</button><div id="description1" class="collapse in" style="position:relative;right:0px;font-size:14px;">This exploratory analysis allows individual variables to be investigated within each dataset, comparing data derived from the HD transgenic (OVT73) sheep to controls. The individual measurements for each animal are visualised using box-plots, with each animals\' unique identifier number displayed. Student\'s t-test is then applied as a measure for comparison (transgenic versus control), producing a p-value to assess the level of significance. In addition, a second p-value for significance is generated, using a linear model to fit sex into the transgenic versus control comparison. A p-value < 0.05 is considered statistically significant at the nominal level.<br><br><div style="font-weight:bold;">Please note: Due to small sample sizes this p-value is exploratory and should only be used as an initial indicator of transgenic versus control differences. It is not adjusted for multiple testing. The next page (Bootstrap and Permutation Tests) can be used to generate a more reliable p-value.</div><br><br>The plots can be downloaded as a PDF by selecting the "Download" button at the bottom of the selection box.<br></div></div>'),
                                      textOutput("Gene_T_test_DL_CT"),
                                      br(),
                                      plotOutput("Gene_Boxplot_CT"),
                                      br(),
                                      br()
                                      #textOutput("Gene_T_test_DL_EvR"),
                                      #br(),
                                      #plotOutput("Gene_Boxplot_ER"),
                                      #br(),
                                      #br()
                                    ))
                         ),
                         # Third tab
                         tabPanel(("Bootstrap and Permutation Tests"),
                                  HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Bootstrap and Permutation Tests</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  sidebarLayout(
                                    sidebarPanel(#Creates a selection box.
                                                 selectInput("select_17", label = h5("Select dataset"),
                                                             choices = list("Transcript expression (nanoString)" = "nano_24", 
                                                                            #"Transcript expression (QPCR)" = "qpcr_all",
                                                                            "Transcript expression (RNAseq Striatum)" = "paul_rnaseq_pval",
                                                                            "Transcript expression (RNAseq Striatal Neurons)" = "neuron_rnaseq",
                                                                            "Proteomics" = "Proteomics",
                                                                            "Metabolites (GC-MS)" = "metab_all",
                                                                            "Metabolites (LC-MS)" = "biocrates_all",
                                                                            "Urea Quantification" = "urea_all")
                                                 ),
                                                 #Creates a selection box.
                                                 selectInput("select_18", label = h5("Select tissue"), 
                                                             choices = list("Dorsolateral" = "nano_24_DL_bs.df", 
                                                                            "Dorsomedial" = "nano_24_DM_bs.df")
                                                 ),
                                                 #Creates a selection box. This takes the names from the dataset columns as the selection variables and the data associated with the variable.
                                                 selectInput("select_19", label = h5("Select specific variable of interest (transcript / metabolite / protein)"), 
                                                             choices = names(nano_24)[6:31]),
                                                 downloadButton("pdflink5")
                                    ),
                                    mainPanel(
                                      HTML('<div class="container2"; style="display:block;width:100%;"><button {class="more/less2" data-toggle="collapse" data-target="#description2" style="position: relative; left:0px; border-radius: 5px; background-color: #00467f;border: none; color: white; padding: 5px 8px; text-align: center; text-decoration: none;font-size: 12px; margin-bottom: 10px; cursor: pointer;"  >More/Less</button><div id="description2" class="collapse in" style="position:relative;right:0px;font-size:14px;">Bootstrap and permutation tests are useful statistical tools that can be performed when sample sizes in the discrete groups are small, thus overcoming the problem of small sample number. Each statistical test can be applied to individual variables of interest within each dataset.<br><br>The bootstrap test is used to identify the probability of observing the true group mean for a variable (the sum of the observed values, divided by the number of samples), given the expression values observed in the sheep. This is done for both the transgenic and control samples. The test relies on random sampling with replacement, and allows the estimation of the sampling distribution. By bootstrapping the expression values and taking the mean 1000 times, a spread of possible group means is presented as a histogram for transgenic (red) and control (blue) groups. The true group mean within each group is displayed on the graph as the solid vertical lines.<br><br>The permutation test calculates a more reliable p-value for expression level comparisons between control and transgenic sheep for the selected variable of interest. The expression values obtained from both the control and transgenic samples are randomly separated out into two separate groups. The difference in expression level between these two groups is then calculated, and repeated 1000 times. The p-value is given by the number of observed differences at least as extreme as the original difference between the control and transgenic sheep, thus giving the significance of the original difference observed.<br><br><div style="font-weight:bold;"> Please note: as the permutation test is run 1000 times, the significance of the p-value is to 0.001.</div><div style="font-weight:bold"></div></div></div>'),
                                      br(),
                                      plotOutput("Bootstrap_plot"),
                                      br(),
                                      br(),
                                      plotOutput("Permutation_plot"),
                                      br(),
                                      br()
                                    )
                                  )),
                         # Fourth tab
                         tabPanel(("PCA"),
                                  HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Principal Components Analysis</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  sidebarLayout(sidebarPanel(# Creates a selection box that allows for multiple selected input variables
                                                             selectInput("select_24", label = h5("Select dataset"), 
                                                                         choices = list("Transcript expression (nanoString) Dorsolateral" = "Dorsolateral", 
                                                                                        "Transcript expression (nanoString) Dorsomedial" = "Dorsomedial",
                                                                                        "Metabolites (GC-MS) Cerebellum",
                                                                                        "Metabolites (GC-MS) Motor Cortex",
                                                                                        "Metabolites (GC-MS) Liver",
                                                                                        "Metabolites (GC-MS) Hippocampus",
                                                                                        "Metabolites (LC-MS) Cerebellum",
                                                                                        "Metabolites (LC-MS) Motor Cortex",
                                                                                        "Metabolites (LC-MS) Liver",
                                                                                        "Metabolites (LC-MS) Plasma",
                                                                                        #"qPCR Motor Cortex",
                                                                                        #"qPCR Anti-striatum",
                                                                                        #"qPCR Foot Plate",
                                                                                        "Urea Serum",
                                                                                        "Urea Urine",
                                                                                        "Urea Cerebellum",
                                                                                        "Urea Hippocampus",
                                                                                        "Urea Motor Cortex",
                                                                                        "Urea Striatum",
                                                                                        "Urea Bladder",
                                                                                        "Urea Heart",
                                                                                        "Urea Kidney",
                                                                                        "Urea Liver",
                                                                                        "Urea Testes",
                                                                                        "Cerebellum Proteomics",
                                                                                        "Motor Cortex Proteomics",
                                                                                        "Striatum Proteomics",
                                                                                        "Striatum RNASeq data",
                                                                                        "Striatal Neuron RNASeq data"
                                                                         ),
                                                                         multiple = TRUE, selected = "Dorsolateral"
                                                             ),
                                                             # Creates a checkbox select input allowing multiple variables to be selected
                                                             checkboxGroupInput("select_25", label = h5("Select which variable arrows to be represented on the 2D PCA plot"), choices = unique(colnames(nano_24_DL)[5:30]), selected = "AP2S1", inline = TRUE),
                                                             # Creates a radio button selection variable where only one selection can be made at a time
                                                             downloadButton("pdflink6")
                                  ),
                                  mainPanel(
                                    HTML('<div class="container3"; style="display:block;width:100%;"><button {class="more/less3" data-toggle="collapse" data-target="#description3" style="position: relative; left:0px; border-radius: 5px; background-color: #00467f;border: none; color: white; padding: 5px 8px; text-align: center; text-decoration: none;font-size: 12px; margin-bottom: 10px; cursor: pointer;"  >More/Less</button><div id="description3" class="collapse in" style="position:relative;right:0px;font-size:14px;">Principal components analysis is one of the most common multivariate dataset tests conducted. This test is used to visualise individual samples in fewer dimensions.<br><br>The first plot on this page shows the level of variance explained between the samples by each of the principal components, (the first 10 are shown). The values for these are given in the text output below this plot. More variance explained by the principal components gives a better representation of the relationship between the samples in the reduced dimensions (given by the cumulative proportion value).<br><br>The relationship between the samples in two-dimensions can be observed in the "Principal Components Plot". Here, the first two principal components are used as the x- and y-axis respectively, and the samples are plotted in the two dimensions, based on levels of variables observed in the datasets used to conduct the test. This visualises the relationship between data derived from the control and HD transgenic sheep.<br><br>This test can be conducted on any number and any combination of datasets as desired by adding to those present in the dataset selection box.<br><br>Another feature present is the availability to observe the effect of each variable on an individual samples position on the two-dimensional plot by selecting the check-box in the selection panel. The arrow added to the plot shows the direction in which an increase in the level of that particular variable will move the sample in the two dimensions. The magnitude of effect of the variable is given by the absolute length of the arrow, thus longer arrows represent a more significant effect on the position of the sample in the two dimensions.<br><br>By default, the Principal Components test is conducted with the data being centred around zero. This is achieved by subtracting the total mean (the mean of all twelve samples, transgenics and controls) for each variable, off the observed value for the variable in each sample.<br><br>Finally, the test is conducted by scaling the variance to a unit value, enabling a more informative comparison between the effects of different variables. This can be de-selected by clicking "No" for the option of whether to scale the variance or not. This is not recommended due to the loss of comparison information at the individual variable level.</div></div>'),
                                    br(),
                                    plotOutput("Scree_plot"),
                                    br(),
                                    br(),
                                    verbatimTextOutput("PCA_summary_stats"),
                                    br(),
                                    br(),
                                    helpText(h4("Top five influential variables for each of the twelve Principal Components")),
                                    verbatimTextOutput("PCA_influential_variables"),
                                    br(),
                                    br(),
                                    plotOutput("first_PCA_Plot"),
                                    br(),
                                    br()
                                  ))
                         ),
                         
                         # Fifth tab
                         tabPanel(("Differential Correlation Statistics"),
                                  HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Differential Correlation Statistics</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  sidebarLayout(
                                    sidebarPanel(selectInput("select_13", label = h5("Select dataset"),
                                                             choices = list("Transcript expression (nanoString)" = "nano_24", 
                                                                            #"Transcript expression (QPCR)" = "qpcr_all",
                                                                            "Transcript expression (RNAseq Striatum)" = "paul_rnaseq_pval",
                                                                            "Transcript expression (RNAseq Striatal Neurons)" = "neuron_rnaseq",
                                                                            "Proteomics" = "Proteomics",
                                                                            "Metabolites (GC-MS)" = "metab_all",
                                                                            "Metabolites (LC-MS)" = "biocrates_all"
                                                                            #"Urea Quantification" = "urea_all"
                                                                            )),
                                                 selectInput("select_14", label = h5("Select tissue"),
                                                             choices = list("Dorsolateral" = "DL_C_and_T_cors_compcor", 
                                                                            "Dorsomedial" = "DM_C_and_T_cors_compcor")),
                                                 downloadButton("text_file")
                                    ),
                                    mainPanel(
                                      HTML('<div class="container7"; style="display:block;width:100%;"><button {class="more/less7" data-toggle="collapse" data-target="#description7" style="position: relative; left:0px; border-radius: 5px; background-color: #00467f;border: none; color: white; padding: 5px 8px; text-align: center; text-decoration: none;font-size: 12px; margin-bottom: 10px; cursor: pointer;"  >More/Less</button><div id="description7" class="collapse in" style="position:relative;right:0px;font-size:14px;">Changes in correlation structure can provide insight into underlying regulatory networks and how these networks may be differentially implicated in a disease process. This analysis presents an exploratory method of investigating the correlation structures between two variables within each dataset, with comparison between transgenic and control groups.<br><br>The analysis determines the Pearson\'s correlation coefficient and associated p-value between two variables (i.e. two genes) for each of the control and transgenic sample groups. It then uses the fisher r-to-z statistic to assess the significance of the difference between the two correlation coefficients.<br><br>This table is ordered by the fisher r-to-z statistic. Only variable-variable combinations with significant fisher-r-to-z statistics (p < 0.05) are displayed. Where applicable, the top fifty correlations are presented, however this may be fewer for some datasets, where the number of significant variable correlations does not reach fifty.<br><br>The search bar above the table enables users to search for specific variables presented in the table.<br><br>This list is downloadable as a .CSV file.</div></div>'),
                                      br(),
                                      DT::dataTableOutput("Top_10_list", width = 300),
                                      br(),
                                      br(),
                                      br()
                                    ))),
                         
                         # Sixth tab
                         tabPanel(("Differential Correlation Plots"),
                                  #Creating the Differential Correlation Plots page.
                                  HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Differential Correlation Plots</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  sidebarLayout(
                                    sidebarPanel(#Creates a selection box for the different datasets.
                                                 selectInput("select_4", label = h5("Select dataset"),
                                                             choices = list("Transcript expression (nanoString)" = "nano_24", 
                                                                            #"Transcript expression (QPCR)" = "qpcr_all",
                                                                            "Transcript expression (RNAseq Striatum)" = "paul_rnaseq_pval",
                                                                            "Transcript expression (RNAseq Striatal Neurons)" = "neuron_rnaseq",
                                                                            "Proteomics" = "Proteomics",
                                                                            "Metabolites (GC-MS)" = "metab_all",
                                                                            "Metabolites (LC-MS)" = "biocrates_all",
                                                                            "Urea Quantification" = "urea_all")),
                                                 #Creates a selection box for the different tissues present in the dataset.
                                                 selectInput("select_5", label = h5("Select tissue"),
                                                             choices = list("Dorsolateral" = "nano_24_DL", 
                                                                            "Dorsomedial" = "nano_24_DM")),
                                                 #Creates a selection box. This takes the names from the dataset columns as the selection variables and the data
                                                 #associated with the variable.
                                                 selectInput("select_6", label = h5("Select variable 1"), 
                                                             choices = names(nano_24)[6:31]),
                                                 #Creates a selection box. This takes the names from the dataset columns as the selection variables and the data
                                                 #associated with the variable.
                                                 selectInput("select_7", label = h5("Select variable 2"),
                                                             choices = names(nano_24)[6:31]),
                                                 downloadButton("pdflink2")
                                    ),
                                    # Show the layout of the page.
                                    mainPanel(
                                      HTML('<div class="container8"; style="display:block;width:100%;"><button {class="more/less8" data-toggle="collapse" data-target="#description8" style="position: relative; left:0px; border-radius: 5px; background-color: #00467f;border: none; color: white; padding: 5px 8px; text-align: center; text-decoration: none;font-size: 12px; margin-bottom: 10px; cursor: pointer;"  >More/Less</button><div id="description8" class="collapse in" style="position:relative;right:0px;font-size:14px;">The differential correlations investigated within each dataset (as explored in the previous tab) can be visualised here. This is another exploratory tool to investigate relationships between different variables, within the control and transgenic groups. Gains or losses in correlation structure may be indicative of differential regulatory processes or associations that occur in the transgenic model. <br><br>Differential correlation plots allow us to visualise variable-variable associations within each group. Each dot represents an individual sample. The correlation coefficient (r) and associated p-value (p) is provided within each plot. The line of best fit is shown in blue for control plots and red for transgenic plots, with confidence intervals in grey. This method of visualisation allows the identification of outliers that cannot be observed in the differential correlation statistics page.<br><br>NOTE: Axes scales differ for each plot.<br><br></div></div>'),
                                      h5(textOutput("CT_comparison")),
                                      verbatimTextOutput("two_variable_C_and_T_comparison"),
                                      br(),
                                      plotOutput("CT_combined_dataset_compcor_plot"),
                                      br(),
                                      #textOutput("ER_comparison"),
                                      #br(),
                                      #verbatimTextOutput("two_variable_E_and_R_comparison"),
                                      #br(),
                                      #plotOutput("ER_combined_dataset_compcor_plot"),
                                      br()
                                    ))),
                         
                         # Seventh tab
                         #tabPanel(("Gene Set Enrichment"),
                                  #Creating the Differential Correlation Plots page.
                        #          HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>Gene Set Enrichment</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                  #Creates a selection box for the different datasets.
                        #          sidebarLayout(
                        #            sidebarPanel(selectInput("select_9", label = h5("Select gene set"),
                      #                                       choices = list("Transcript expression (nanoString)" = "nanoString_DAVID_analysis",
                      #                                                      "Transcript expression (RNAseq Striatum)" = "Striatum_RNASeq_DAVID_analysis",
                      #                                                      "Transcript expression (RNAseq Striatal Neurons)" = "Striatal_Neuron_DAVID_analysis")),
                      #                           downloadButton("pdflink7")),
                                    # Show the layout of the page.
                      #              mainPanel(
                      #                HTML('<div class="container9"; style="display:block;width:100%;right:0px;"><button {class="more/less9" data-toggle="collapse" data-target="#description9" style="position: relative; right:0px; border-radius: 5px; background-color: #00467f;border: none; color: white; padding: 5px 8px; text-align: center; text-decoration: none;font-size: 12px; margin-bottom: 10px; cursor: pointer;"  >More/Less</button><div id="description9" class="collapse in" style="position:relative;right:0px;font-size:14px;">Gene Set Enrichment analysis is a method to identify classes of genes that are overrepresented in a larger set of genes.<br><br>Three transcriptomic datasets (targeted gene expression by nanoString, Striatum RNASeq, and laser captured Striatal Neuron RNASeq) initially underwent differential expression analyses, yielding lists of genes that were significantly up- or downregulated in the OVT73 animals compared to controls (p<0.05). These gene lists were then input into the Database for Annotation, Visualisation and Integrated Discovery (DAVID) Bioinformatics Resource 6.8 for gene enrichment analysis. Annotation categories were pre-selected for each enrichment analysis. The default gene ontology categories were \'GOTERM_BP_DIRECT\', \'GOTERM_CC_DIRECT\' and \'GOTERM_MF_DIRECT\' for biological process, cellular component and molecular function respectively, along with one pathway category; \'KEGG_PATHWAY\'. Since Huntington\'s disease is a human disease, annotations were limited to those identified in \'Homo sapiens\' under the species categorisation.<br><br>The output from the DAVID analysis for each of the datasets is presented in a table. The term, pathway and implicated genes are given, along with statistical values to determine the strength of the relationship of the genes within the pathway. The search bar above the table enables users to search for specific pathways or genes presented in the table.</div></div>'),
                      #                br(),
                      #                DT::dataTableOutput("Gene_set_table", width = 300),
                      #                br(),
                      #                br(),
                      #                br()
                      #                )
                      #            )),
                         
                         
                       # Eighth tab - create a second drop-down select tab in the navbar
                                  tabPanel(("About/Contact"),
                                           HTML("<div class='wrap' style='width:100%;height:111px;overflow:hidden;background-color:white;'><h1 class='header' style='margin-top:20px;position:absolute;display:inline-block;'>About this Project / Contact Us</h1><div class='logo' style='position:relative;overflow: hidden;display:inline-block;float:right;'><img class='icon-overlay' src='HD-brain-overlay_sheep_2.png' style='position: relative; width: 96px; height: 96px; z-index: 100;'/><img class='icon-underlay' src='HD-brain-underlay.svg' style='position: absolute; top: 0; left: 0; width: 96px; height: 96px;'/></div></div>"),
                                           sidebarLayout(
                                             sidebarPanel(helpText(h3("Contact us:")),
                                                          textInput("First_Name", label = h5("First Name"), value = ""),
                                                          textInput("Last_Name", label = h5("Last Name"), value = ""),
                                                          textInput("Email", label = h5("Email"), value = ""),
                                                          textAreaInput("Message", label = h5("Comment"), value = "", height = "200px"),
                                                          actionButton("Send", label = "Send")), position = "right",
                                             mainPanel(
                                               tags$head(tags$script(src = "message-handler.js")),
                                               HTML('<div style="display:block;width:100%;"><div style="position:relative;right:0px;font-size:14px;">Huntington\'s Disease (HD) is relatively uncommon, affecting approximately 1 in 10,000 individuals of European origin. The disease is caused by the expansion of a coding polymorphic CAG trinucleotide repeat located in exon 1 of the Huntingtin (HTT) gene <a href = "https://www.omim.org/entry/613004"> OMIM # 613004</a>. The biological functions of HTT, and the pathogenic mechanism mediated by the mutant allele are not yet fully understood. There is no therapy in clinical use that can prevent or delay the onset of HD.<br><br>There have been many animal models of HD made to enable the investigation of the disease process, and also for use in pre-clinical pharmacological testing. A collaborative project between The University of Auckland, Harvard Medical School and the South Australian Research and Development Institute, resulted in the production first large mammalian model of HD - a transgenic sheep line termed OVT73. The OVT73 sheep line carries copies of full length human huntingtin cDNA with an expanded polyglutamine coding repeat of 73 units. The repeat is relatively stable on transmission with an expression level of approximately half an allele. Sheep have many advantages for neurological study and drug testing due to their large size, long lifespan, and complex brain structure which is comparable to humans. The transgenic sheep show no overt neurological symptoms (some >10 years of age), but do develop some of the hallmark brain pathology of HD such as huntingtin positive inclusions and altered expression of genes that are implicated in HD. The sheep have a measurable circadian abnormality, which is a characteristic behaviour of HD patients and a metabolic disruption. There are approximately one thousand OVT73 and control sheep available for research on the farm which is a normal pasture-based environment. Results from a wide range of OVT73 studies suggest that the model recapitulates the prodromal phase of Huntington\'s disease before motor age of onset.<br><br>A large series of data has been generated from a single cohort of OVT73 animals (n=6) and control sheep (n=6) killed at 5 years of age. This includes transcriptomic, proteomic and metabolic data collected from multiple tissue samples through various independent experimental analyses. This data has been collated and integrated into a single multidimensional platform presented on this website. The aim is to facilitate the \'multi-omic\' exploration of the data for correlations and provides a more holistic view of the OVT73 sheep compared to controls. <br><br>Our aim is to share the OVT73 data with the HD research community, through an interactive database that can be queried for specific genes/proteins/metabolites of interest. The data can be explored using a range of statistical techniques and multivariate methods (presented as separate tabs at the top of the page). Alternatively, the raw datasets can be downloaded for further analyses here. The ultimate goal of this database is to offer another layer of HD data, from sheep, in an attempt to discover the causative mechanisms behind HD progression and identify potential therapeutic targets. <br><br> This interactive, web-based application allows the investigation and exploration of \'multi-omic\' variables expressed in a sheep model of Huntington\'s disease (HD), termed OVT73. OVT73, which has been extensively characterised, represents a prodromal form of Huntington\'s disease (HD). The aim of this application is to share OVT73 data with the wider HD research community; allow variables of interest to be queried in a HD sheep model; bridge the gap between mouse and human HD studies and act as a hypothesis generator for HD mechanisms of pathogenesis.<br><br>If you would like to contact the Snell lab group or provide feedback, please use the Contact Us panel. If you would like to remain anonymous, please leave the name fields blank. This is an early stage project and any feedback provided will be greatly appreciated.<br><br>We would like to thank The University of Auckland for providing the services necessary for hosting and maintaining the application and providing the means to undertake the research that has been conducted.<br><br>We would also like to thank our current and previous collaborators for the opportunity to conduct experiments and allow the publication of resulting data. <br><br>The development of this application was funded by Brain Research New Zealand (BRNZ).<br><br>Created by Matthew Grant and Emily Mears, from the Snell lab group, at The University of Auckland.</div></div>'),
                                               br(),
                                               br())
                                           ),
                                           HTML('<div class="wrap" style="width:100%;height:150px;background-color:white;"> <img src="BRNZ.png" style="width:100px;height:100px;display:inline-block;margin-top:25px;margin-left:15px"> <img src="UOA.png" style="width:220px;height:100px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="CHDI.png" style="width:180px;height:100px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="SARDI.png" style="width:60px;height:100px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="UniofAdelaide.png" style="width:160px;height:100px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="SAHMRI.png" style="width:180px;height:100px;display:inline-block;margin-top:25px;margin-left:1%"> <img src="CBR.png" style="width:220px;height:100px;display:inline-block;margin-top:25px;margin-left:1%"> </div>')
                                  )
    # Tab 8B
                                    
                         )
)

# This is the end of the User Interface initialisation
# The following code is that which is used to create the graphics and text that fill the pages in their defined locations,
# which are defined by their name and location in the mainPanel function





##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

# Define server logic required to produce output defined.
server <- shinyServer(function(input, output, session) {c(
  
  ##################################################################################################################################################
  ##################################################################################################################################################
  
  ############ Single variable basic statistics ############
  
  # Creating the drop-down selection box options, which react to the name that is contained in the variable: 'dataset_page_1'
  # These options also relate to a specific dataset
  observe({
    dataset_page_1 <- input$select_1
    if (dataset_page_1 == "metab_all")
      updateSelectInput(session, "select_2",
                        choices = list("Cerebellum" = "metab_cb",
                                       "Motor Cortex" = "metab_mctx",
                                       "Liver" = "metab_liv",
                                       "Hippocampus" = "metab_hipp"))
    if (dataset_page_1 == "nano_24")
      updateSelectInput(session, "select_2",
                        choices = list("Dorsolateral" = "nano_24_DL", "Dorsomedial" = "nano_24_DM"))
    if (dataset_page_1 == "biocrates_all")
      updateSelectInput(session, "select_2",
                        choices = list("Cerebellum" = "biocrates_cb",
                                       "Motor Cortex" = "biocrates_mctx",
                                       "Liver" = "biocrates_liv",
                                       "Plasma" = "biocrates_plasma"))
    if (dataset_page_1 == "qpcr_all")
      updateSelectInput(session, "select_2",
                        choices = list("Motor Cortex" = "qpcr_mctx",
                                       "Anti-striatum" = "qpcr_antistr",
                                       "Floor Plate" = "qpcr_fp"))
    if (dataset_page_1 == "urea_all")
      updateSelectInput(session, "select_2",
                        choices = list("Serum" = "urea_serum",
                                       "Urine" = "urea_urine",
                                       "Cerebellum" = "urea_cb",
                                       "Hippocampus" = "urea_hipp",
                                       "Motor Cortex" = "urea_mctx",
                                       "Striatum" = "urea_stri",
                                       "Bladder" = "urea_blad",
                                       "Heart" = "urea_heart",
                                       "Kidney" = "urea_kid",
                                       "Liver" = "urea_liv",
                                       "Testes" = "urea_testes"))
    if (dataset_page_1 == "Proteomics")
      updateSelectInput(session, "select_2",
                        choices = list("Cerebellum" = "proteo_cb_sig_only",
                                       "Motor Cortex" = "proteo_mctx_sig_only",
                                       "Striatum" = "proteo_str_sig_only"))
    if (dataset_page_1 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_2",
                        choices = list("Striatum" = "paul_rnaseq_pval"))
    if (dataset_page_1 == "neuron_rnaseq")
      updateSelectInput(session, "select_2",
                        choices = list("Striatal Neurons" = "neuron_rnaseq"))
    
    
    
  }),
  # Creating the names that appear in the drop-down selection box to select the variable names present in the selected dataset from the previous selection box
  observe({
    dataset_page_1 <- input$select_1
    if (dataset_page_1 == "metab_all")
      updateSelectInput(session, "select_3",
                        choices = names(metab_all)[10:71])
    if (dataset_page_1 == "nano_24")
      updateSelectInput(session, "select_3",
                        choices = names(nano_24)[6:31])
    if (dataset_page_1 == "biocrates_all")
      updateSelectInput(session, "select_3",
                        choices = names(biocrates_all[6:173]))
    if (dataset_page_1 == "qpcr_all")
      updateSelectInput(session, "select_3",
                        choices = names(qpcr_all[10:19]))
    if (dataset_page_1 == "urea_all")
      updateSelectInput(session, "select_3",
                        choices = names(urea_all[8:9]))
    if (input$select_2 == "proteo_cb_sig_only")
      updateSelectInput(session, "select_3",
                        choices = names(proteo_cb_sig_only[5:25]))
    if (input$select_2 == "proteo_mctx_sig_only")
      updateSelectInput(session, "select_3",
                        choices = names(proteo_mctx_sig_only[5:21]))
    if (input$select_2 == "proteo_str_sig_only")
      updateSelectInput(session, "select_3",
                        choices = names(proteo_str_sig_only[5:43]))
    if (dataset_page_1 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_3",
                        choices = names(paul_rnaseq_pval[3:819]))
    if (dataset_page_1 == "neuron_rnaseq")
      updateSelectInput(session, "select_3",
                        choices = names(neuron_rnaseq[3:47]))
    
  }),
  
  # Retrieving the dataset that is specified by the selection of 'input$select_2'
  page_1_df <- reactive({
    w <- get(input$select_2)
  }),
  # Creating the plot that will be positioned in the 'Gene_Boxplot_CT' section of the main panel of tab 2
  output$Gene_Boxplot_CT <- renderPlot({
    # Retrieve the reactive variable and store in a usable variable
    p1.df = page_1_df()
    # Creating the plot within the renderPlot function
    # Using the status of the sheep to separate the individuals into two groups (x-axis) and selecting out the specific variable that is to be investigated, specified by 'input$select_3' (y-axis)
    # Other information is added using the functions of ggplot2 including the titles, colours, and names.
    ggplot(p1.df, aes(x = p1.df$status, y = p1.df[,which(colnames(p1.df) == input$select_3)], label = rownames(p1.df), colour = p1.df$status, fill = p1.df$status)) + geom_boxplot(show.legend = FALSE, colour = "black", alpha=I(.5)) + ggtitle(paste(input$select_3, "expression in", p1.df$tissue[1], "by transgene status")) + labs(subtitle = paste("T-test p-value for", input$select_3, "by transgene status:", round(t.test(p1.df[p1.df$status == "C", colnames(p1.df) == input$select_3], p1.df[p1.df$status == "T", colnames(p1.df) == input$select_3], var.equal = TRUE)$p.value,4), "\nLinear model p-value significance of sex effect:", round(summary(lm(p1.df[,which(colnames(p1.df) == input$select_3)] ~ p1.df$status * p1.df$sex))[[4]][15], 4))) + theme(plot.title=element_text(size=15, face="bold")) + theme(plot.subtitle=element_text(size=12, face="italic")) + geom_text(aes(colour = factor(p1.df$status)), show.legend = FALSE, fontface = "bold.italic", size = 5) + ylab(paste(input$select_3, "level")) + xlab("Status") + scale_x_discrete(labels = c("Control", "Transgenic")) + theme_bw() + scale_fill_manual(values = c("blue3", "red3")) + scale_colour_manual(values = c("red3", "blue3")) + theme(text = element_text(size = 15, face = "italic"))
  }),
  
  # This code is similar to that above, however the split of the individuals is by sex - to investigate potential sex dependent traits
  #output$Gene_Boxplot_ER <- renderPlot({
  #  p1.df = page_1_df()
  #  ggplot(p1.df, aes(x = p1.df$sex, y = p1.df[,which(colnames(p1.df) == input$select_3)], label = rownames(p1.df), colour = p1.df$sex, fill = p1.df$sex)) + geom_boxplot(show.legend = FALSE, colour = "black", alpha=I(.5)) + ggtitle(paste(input$select_3, "expression in", p1.df$tissue[1], "by sex")) + labs(subtitle = paste("T-test p-value for", input$select_3, "by sex:", round(t.test(p1.df[p1.df$sex == "E", colnames(p1.df) == input$select_3], p1.df[p1.df$sex == "R", colnames(p1.df) == input$select_3], var.equal = TRUE)$p.value,4))) + theme(plot.title=element_text(size=15, face="bold")) + theme(plot.subtitle=element_text(size=12, face="italic")) + geom_text(aes(colour = factor(p1.df$sex)), show.legend = FALSE, fontface = "bold.italic", size = 5) + ylab(paste(input$select_3, "expression")) + xlab("Sex") + scale_x_discrete(labels = c("Ewe", "Ram")) + theme_bw() + scale_fill_manual(values = c("blue3", "red3")) + scale_colour_manual(values = c("red3", "blue3")) + theme(text = element_text(size = 15, face = "italic"))
  #}),
  
  # The plots presented on the page are downloadable - this code creates the plots in a PDF format that is downloaded when the button is selected
  # First the name of the file is created, this includes the information of tissue type and name of the variable selected
  output$pdflink1 <- downloadHandler(
    filename = function(){
      page_1_df <- reactive({
        w <- get(input$select_2)
      })
      p1.df = page_1_df()
      paste(p1.df$tissue[[1]], input$select_3, "boxplot.pdf")
    },
    # Below, the two plots are created again in the PDF environment
    content = function(file) {
      page_1_df <- reactive({
        w <- get(input$select_2)
      })
      p1.df = page_1_df()
      pdf(file, width = 10, height = 10)
      p1.df = page_1_df()
      grid.arrange(
        ggplot(p1.df, aes(x = p1.df$status, y = p1.df[,which(colnames(p1.df) == input$select_3)], label = rownames(p1.df), colour = p1.df$status, fill = p1.df$status)) + geom_boxplot(show.legend = FALSE, colour = "black", alpha=I(.5)) + ggtitle(paste(input$select_3, "expression in", p1.df$tissue[1], "by transgene status")) + labs(subtitle = paste("T-test p-value for", input$select_3, "by transgene status:", round(t.test(p1.df[p1.df$status == "C", colnames(p1.df) == input$select_3], p1.df[p1.df$status == "T", colnames(p1.df) == input$select_3], var.equal = TRUE)$p.value,4), "\nLinear Model p-value for sex effect:", round(summary(lm(p1.df[,which(colnames(p1.df) == input$select_3)] ~ p1.df$status * p1.df$sex))[[4]][15], 4))) + theme(plot.title=element_text(size=15, face="bold")) + theme(plot.subtitle=element_text(size=12, face="italic")) + geom_text(aes(colour = factor(p1.df$status)), show.legend = FALSE, fontface = "bold.italic", size = 5) + ylab(paste(input$select_3, "expression")) + xlab("Status") + scale_x_discrete(labels = c("Control", "Transgenic")) + theme_bw() + scale_fill_manual(values = c("blue3", "red3")) + scale_colour_manual(values = c("red3", "blue3")) + theme(text = element_text(size = 15, face = "italic")),
        ncol = 1, nrow = 1)
      dev.off()
    }
  ),
  
  ###############################################################################################################################################################
  ########################################################### Bootstrap and Permutation Tests ###################################################################
  ###############################################################################################################################################################
  
  # Creating the drop down selection boxes
  observe({
    dataset_page_6 <- input$select_17
    if (dataset_page_6 == "metab_all")
      updateSelectInput(session, "select_18",
                        choices = list("Cerebellum" = "metab_cb_bs.df",
                                       "Motor Cortex" = "metab_mctx_bs.df",
                                       "Liver" = "metab_liv_bs.df",
                                       "Hippocampus" = "metab_hipp_bs.df"))
    if (dataset_page_6 == "nano_24")
      updateSelectInput(session, "select_18",
                        choices = list("Dorsolateral" = "nano_24_DL_bs.df",
                                       "Dorsomedial" = "nano_24_DM_bs.df"))
    if (dataset_page_6 == "biocrates_all")
      updateSelectInput(session, "select_18",
                        choices = list("Cerebellum" = "biocrates_cb_bs.df",
                                       "Motor Cortex" = "biocrates_mctx_bs.df",
                                       "Liver" = "biocrates_liv_bs.df",
                                       "Plasma" = "biocrates_plasma_bs.df"))
    if (dataset_page_6 == "qpcr_all")
      updateSelectInput(session, "select_18",
                        choices = list("Motor Cortex" = "qpcr_mctx_bs.df",
                                       "Anti-striatum" = "qpcr_antistr_bs.df",
                                       "Floor Plate" = "qpcr_fp_bs.df"))
    if (dataset_page_6 == "urea_all")
      updateSelectInput(session, "select_18",
                        choices = list("Serum" = "urea_serum_bs.df",
                                       "Urine" = "urea_urine_bs.df",
                                       "Cerebellum" = "urea_cb_bs.df",
                                       "Hippocampus" = "urea_hipp_bs.df",
                                       "Motor Cortex" = "urea_mctx_bs.df",
                                       "Striatum" = "urea_stri_bs.df",
                                       "Bladder" = "urea_blad_bs.df",
                                       "Heart" = "urea_heart_bs.df",
                                       "Kidney" = "urea_kid_bs.df",
                                       "Liver" = "urea_liv_bs.df",
                                       "Testes" = "urea_testes_bs.df"))
    if (dataset_page_6 == "Proteomics")
      updateSelectInput(session, "select_18",
                        choices = list("Cerebellum" = "proteo_cb_sig_only_bs.df",
                                       "Motor Cortex" = "proteo_mctx_sig_only_bs.df",
                                       "Striatum" = "proteo_str_sig_only_bs.df"))
    if (dataset_page_6 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_18",
                        choices = list("Striatum" = "paul_rnaseq_pval_bs.df"))
    if (dataset_page_6 == "neuron_rnaseq")
      updateSelectInput(session, "select_18",
                        choices = list("Striatal Neurons" = "neuron_rnaseq_bs.df"))
    
  }),
  
  observe({
    dataset_page_6 <- input$select_17
    if (dataset_page_6 == "metab_all")
      updateSelectInput(session, "select_19",
                        choices = names(metab_all)[10:71])
    if (dataset_page_6 == "nano_24")
      updateSelectInput(session, "select_19",
                        choices = names(nano_24)[6:31])
    if (dataset_page_6 == "biocrates_all")
      updateSelectInput(session, "select_19",
                        choices = names(biocrates_all[6:173]))
    if (dataset_page_6 == "qpcr_all")
      updateSelectInput(session, "select_19",
                        choices = names(qpcr_all[10:19]))
    if (dataset_page_6 == "urea_all")
      updateSelectInput(session, "select_19",
                        choices = names(urea_all[8:9]))
    if (input$select_18 == "proteo_cb_sig_only_bs.df")
      updateSelectInput(session, "select_19",
                        choices = names(proteo_cb_sig_only[5:25]))
    if (input$select_18 == "proteo_mctx_sig_only_bs.df")
      updateSelectInput(session, "select_19",
                        choices = names(proteo_mctx_sig_only[5:21]))
    if (input$select_18 == "proteo_str_sig_only_bs.df")
      updateSelectInput(session, "select_19",
                        choices = names(proteo_str_sig_only[5:43]))
    if (dataset_page_6 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_19",
                        choices = names(paul_rnaseq_pval[3:819]))
    if (dataset_page_6 == "neuron_rnaseq")
      updateSelectInput(session, "select_19",
                        choices = names(neuron_rnaseq[3:47]))
    
  }),
  
  # Retrieving the relevant datasets
  page_6_df <- reactive({
    w <- get(input$select_18)
  }),
  
  page_6.2_df <- reactive({
    z <- get(input$select_17)
  }),
  
  output$Bootstrap_plot <- renderPlot({
    p6.df = page_6_df()
    ggplot(p6.df, aes_string(x=input$select_19)) +
      geom_density(fill = "Blue", alpha = 0.4) +
      geom_density(aes_string(x=paste(input$select_19, ".1", sep = "")),fill = "Red", alpha = 0.4) + 
      theme_bw() + 
      geom_vline(aes(xintercept = (mean(data_sets_list2[[input$select_18]][7:12, colnames(data_sets_list2[[input$select_18]]) == input$select_19]))), col = "Red", size = 1) +
      geom_vline(aes(xintercept = (mean(data_sets_list2[[input$select_18]][1:6, colnames(data_sets_list2[[input$select_18]]) == input$select_19]))), col = "Blue", size = 1) + 
      labs(title = paste("Bootstrap for control and transgenic", input$select_19, "expression levels"), subtitle = paste("Observed mean for control samples = ", round(mean(data_sets_list2[[input$select_18]][1:6, colnames(data_sets_list2[[input$select_18]]) == input$select_19]),2), "\nObserved mean for transgenic samples = ", round(mean(data_sets_list2[[input$select_18]][7:12, colnames(data_sets_list2[[input$select_18]]) == input$select_19]),2)), x = paste("Expression level of", input$select_19)) + 
      theme(plot.title=element_text(size=18, face="bold"), plot.subtitle=element_text(size=15, face="italic"))
  }),
  
  output$Permutation_plot <- renderPlot({
    qplot(data_sets_list3[[input$select_18]][,colnames(data_sets_list3[[input$select_18]]) == input$select_19],
          geom="histogram",
          main=paste("Permutation for",input$select_19, "expression levels"),
          xlab="Expression level",  
          alpha=I(.5)) + theme_bw() + geom_vline(aes(xintercept=data_sets_list3[[input$select_18]][1000,colnames(data_sets_list3[[input$select_18]]) == input$select_19]), col = "red3", size = 1) + labs(subtitle = paste("Permutation P-value for", input$select_19, "=", 
                                                                                                                                                                                                                           if(data_sets_list3[[input$select_18]][1000, colnames(data_sets_list3[[input$select_18]]) == input$select_19] > 0){
                                                                                                                                                                                                                             ((sum(data_sets_list3[[input$select_18]][1000, colnames(data_sets_list3[[input$select_18]]) == input$select_19] <= data_sets_list3[[input$select_18]][1:999,colnames(data_sets_list3[[input$select_18]]) == input$select_19])+1) / 1000)
                                                                                                                                                                                                                           }
                                                                                                                                                                                                                           else{
                                                                                                                                                                                                                             ((sum(data_sets_list3[[input$select_18]][1000, colnames(data_sets_list3[[input$select_18]]) == input$select_19] >= data_sets_list3[[input$select_18]][1:999,colnames(data_sets_list3[[input$select_18]]) == input$select_19])+1) / 1000)
                                                                                                                                                                                                                           }), size = 5) + theme(plot.title=element_text(size=18, face="bold")) + theme(plot.subtitle=element_text(size=15, face="italic"))
  }),
  
  output$pdflink5 <- downloadHandler(
    filename = function(){
      page_6_df <- reactive({
        w <- get(input$select_18)
      })
      page_6.2_df <- reactive({
        z <- get(input$select_17)
      })
      
      paste(input$select_19, "bootstrap_and_permutations.pdf")
    },
    content = function(file) {
      page_6_df <- reactive({
        w <- get(input$select_18)
      })
      p6.df = page_6_df()
      pdf(file, width = 10, height = 10)
      
      
      grid.arrange(
        ggplot(p6.df, aes_string(x=input$select_19)) +
          geom_density(fill = "blue", alpha = 0.4) +
          geom_density(aes_string(x=paste(input$select_19, ".1", sep = "")),fill = "red", alpha = 0.4) + 
          theme_bw() + 
          geom_vline(aes(xintercept = (mean(data_sets_list2[[input$select_18]][7:12, colnames(data_sets_list2[[input$select_18]]) == input$select_19]))), col = "Red", size = 1) +
          geom_vline(aes(xintercept = (mean(data_sets_list2[[input$select_18]][1:6, colnames(data_sets_list2[[input$select_18]]) == input$select_19]))), col = "Blue", size = 1) + 
          labs(title = paste("Bootstrap for control and transgenic", input$select_19, "expression levels"), subtitle = paste("Observed mean for control samples = ", round(mean(data_sets_list2[[input$select_18]][1:6, colnames(data_sets_list2[[input$select_18]]) == input$select_19]),2), "\nObserved mean for transgenic samples = ", round(mean(data_sets_list2[[input$select_18]][7:12, colnames(data_sets_list2[[input$select_18]]) == input$select_19]),2)), x = paste("Expression level of", input$select_19)) + 
          theme(plot.title=element_text(size=18, face="bold"), plot.subtitle=element_text(size=15, face="italic")),
        
        qplot(data_sets_list3[[input$select_18]][,colnames(data_sets_list3[[input$select_18]]) == input$select_19],
              geom="histogram",
              main=paste("Permutation for",input$select_19, "expression levels"),
              xlab="Expression level", 
              alpha=I(.75)) + geom_vline(aes(xintercept=data_sets_list3[[input$select_18]][1000,colnames(data_sets_list3[[input$select_18]]) == input$select_19]), col = "red3", size = 1) + labs(subtitle = paste("Permutation P-value for", input$select_19, "=", 
                                                                                                                                                                                                                   if(data_sets_list3[[input$select_18]][1000, colnames(data_sets_list3[[input$select_18]]) == input$select_19] > 0){
                                                                                                                                                                                                                     ((sum(data_sets_list3[[input$select_18]][1000, colnames(data_sets_list3[[input$select_18]]) == input$select_19] <= data_sets_list3[[input$select_18]][1:999,colnames(data_sets_list3[[input$select_18]]) == input$select_19])+1) / 1000)
                                                                                                                                                                                                                   }
                                                                                                                                                                                                                   else{
                                                                                                                                                                                                                     ((sum(data_sets_list3[[input$select_18]][1000, colnames(data_sets_list3[[input$select_18]]) == input$select_19] >= data_sets_list3[[input$select_18]][1:999,colnames(data_sets_list3[[input$select_18]]) == input$select_19])+1) / 1000)
                                                                                                                                                                                                                   }), size = 5) + theme(plot.title=element_text(size=15, face="bold")) + theme(plot.subtitle=element_text(size=12, face="italic"))
      )
      dev.off()
    }
  ),
  
  ###############################################################################################################################################################
  ###############################################################################################################################################################
  
  ############ Principal Components Analysis ############
  
  observe({
    updateCheckboxGroupInput(session, "select_25",
                             choices = unique(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))), inline = TRUE)
  }),
  
  output$Scree_plot <- renderPlot({
    screeplot(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE), main = "Scree plot of variance explained by each PC", xlab = "Principal components")
    axis(1, at = seq(0.7, 12, 1.2), tick = TRUE, labels = paste("PC ", 1:10))
  }),
  
  output$PCA_summary_stats <- renderPrint({
    summary(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE))
  }),
  
  output$PCA_influential_variables <- renderPrint({
    cbind("PC1" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,1]))])[1:5],
          "PC2" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,2]))])[1:5],
          "PC3" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,3]))])[1:5],
          "PC4" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,4]))])[1:5],
          "PC5" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,5]))])[1:5],
          "PC6" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,6]))])[1:5],
          "PC7" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,7]))])[1:5],
          "PC8" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,8]))])[1:5],
          "PC9" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,9]))])[1:5],
          "PC10" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,10]))])[1:5],
          "PC11" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,11]))])[1:5],
          "PC12" = colnames(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))))))[,order(abs(prcomp(cbind("status" = factor(c(rep(1, 6), rep(2, 6))), "sex" = factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12, dimnames = list(NULL, as.vector(unlist(sapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names)))))), scale. = TRUE)[[2]][,12]))])[1:5])
  }),
  
  output$first_PCA_Plot <- renderPlot({
    plot(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,1], prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2], pch = NA, xlab = "PC1", ylab = "PC2", main = "Principal Components Plot", bty = 'l')
    abline(h = 0)
    abline(v = 0)
    par(font = 4)
    grid(lwd = 1, lty = "dashed", col = "grey")
    text(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,1], prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2], labels = rownames(nano_24_DL), col = c("Blue4", "Blue4", "Blue4", "Blue4", "Blue4", "Blue4", "Red3", "Red3", "Red3", "Red3", "Red3", "Red3")) 
    abline(h=0,v=0,lty=2,col="grey")
    if(is.null(input$select_25) == FALSE){
      arrows(x0 = 0, y0 = 0,
             x1 = prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,1][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)] * 
               min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))),
             y1 = prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)] * 
               min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))))
      text(x = prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,1][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)] * 
             min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))),
           y = (prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)]) * 
             min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))), labels = input$select_25)
    }
  }),
  
  output$pdflink6 <- downloadHandler(
    filename = function(){
      paste(input$select_24, "PCA plots.pdf")
    },
    content = function(file){
      pdf(file, width = 7.5, height = 5)
      par(mfrow=c(1, 1))
      screeplot(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE), main = "Scree plot of variance explained by each PC", xlab = "Principal components")
      axis(1, at = seq(0.7, 12, 1.2), tick = TRUE, labels = paste("PC ", 1:10))
      plot(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,1], prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2], pch = NA, xlab = "PC1", ylab = "PC2", main = "Principal Components Plot")
      par(font = 4)
      text(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,1], prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2], labels = rownames(nano_24_DL), col = c("Blue4", "Blue4", "Blue4", "Blue4", "Blue4", "Blue4", "Red3", "Red3", "Red3", "Red3", "Red3", "Red3")) 
      abline(h=0,v=0,lty=2,col="grey")
      if(is.null(input$select_25) == FALSE){
        arrows(x0 = 0, y0 = 0,
               x1 = prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,1][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)] * 
                 min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))),
               y1 = prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)] * 
                 min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))))
        text(x = prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,1][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)] * 
               min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))),
             y = (prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2][which(as.vector(unlist(lapply(data_sets_list[which(names(data_sets_list) %in% input$select_24)], names))) %in% input$select_25)]) * 
               min((max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,2])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,2]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]))), (max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[5]][,3])/(max(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3]) - min(prcomp(cbind(factor(c(rep(1, 6), rep(2, 6))), factor(c(1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2)), matrix(unlist(data_sets_list[which(names(data_sets_list) %in% input$select_24)]), nrow = 12)), scale. = TRUE)[[2]][,3])))), labels = input$select_25)
      }
      dev.off()
    }),
  
  
  
  ###############################################################################################################################################################
  ###############################################################################################################################################################
  
  ############ Two variable correlations ############
  
  observe({
    dataset_page_2 <- input$select_4
    if (dataset_page_2 == "metab_all")
      updateSelectInput(session, "select_5",
                        choices = list("Cerebellum" = "metab_cb",
                                       "Motor Cortex" = "metab_mctx",
                                       "Liver" = "metab_liv",
                                       "Hippocampus" = "metab_hipp"))
    if (dataset_page_2 == "nano_24")
      updateSelectInput(session, "select_5",
                        choices = list("Dorsolateral" = "nano_24_DL",
                                       "Dorsomedial" = "nano_24_DM"))
    
    if (dataset_page_2 == "biocrates_all")
      updateSelectInput(session, "select_5",
                        choices = list("Cerebellum" = "biocrates_cb",
                                       "Motor Cortex" = "biocrates_mctx",
                                       "Liver" = "biocrates_liv",
                                       "Plasma" = "biocrates_plasma"))
    
    if (dataset_page_2 == "qpcr_all")
      updateSelectInput(session, "select_5",
                        choices = list("Motor Cortex" = "qpcr_mctx",
                                       "Anti-striatum" = "qpcr_antistr",
                                       "Floor Plate" = "qpcr_fp"))
    
    if (dataset_page_2 == "urea_all")
      updateSelectInput(session, "select_5",
                        choices = list("Serum" = "urea_serum",
                                       "Urine" = "urea_urine",
                                       "Cerebellum" = "urea_cb",
                                       "Hippocampus" = "urea_hipp",
                                       "Motor Cortex" = "urea_mctx",
                                       "Striatum" = "urea_stri",
                                       "Bladder" = "urea_blad",
                                       "Heart" = "urea_heart",
                                       "Kidney" = "urea_kid",
                                       "Liver" = "urea_liv",
                                       "Testes" = "urea_testes"))
    if (dataset_page_2 == "Proteomics")
      updateSelectInput(session, "select_5",
                        choices = list("Cerebellum" = "proteo_cb_sig_only",
                                       "Motor Cortex" = "proteo_mctx_sig_only",
                                       "Striatum" = "proteo_str_sig_only"))
    if (dataset_page_2 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_5",
                        choices = list("Striatum" = "paul_rnaseq_pval"))
    if (dataset_page_2 == "neuron_rnaseq")
      updateSelectInput(session, "select_5",
                        choices = list("Striatal Neurons" = "neuron_rnaseq"))
    
  }),
  
  observe({
    dataset_page_2 <- input$select_4
    if (dataset_page_2 == "metab_all")
      updateSelectInput(session, "select_6",
                        choices = names(metab_all)[10:71])
    if (dataset_page_2 == "nano_24")
      updateSelectInput(session, "select_6",
                        choices = names(nano_24)[6:31])
    if (dataset_page_2 == "biocrates_all")
      updateSelectInput(session, "select_6",
                        choices = names(biocrates_all[6:173]))
    if (dataset_page_2 == "qpcr_all")
      updateSelectInput(session, "select_6",
                        choices = names(qpcr_all[10:19]))
    if (dataset_page_2 == "urea_all")
      updateSelectInput(session, "select_6",
                        choices = names(urea_all[8:9]))
    if (input$select_5 == "proteo_cb_sig_only")
      updateSelectInput(session, "select_6",
                        choices = names(proteo_cb_sig_only[5:25]))
    if (input$select_5 == "proteo_mctx_sig_only")
      updateSelectInput(session, "select_6",
                        choices = names(proteo_mctx_sig_only[5:21]))
    if (input$select_5 == "proteo_str_sig_only")
      updateSelectInput(session, "select_6",
                        choices = names(proteo_str_sig_only[5:43]))
    if (dataset_page_2 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_6",
                        choices = names(paul_rnaseq_pval[3:819]))
    if (dataset_page_2 == "neuron_rnaseq")
      updateSelectInput(session, "select_6",
                        choices = names(neuron_rnaseq[3:47]))
    
    
  }),
  
  observe({
    dataset_page_2 <- input$select_4
    if (dataset_page_2 == "metab_all")
      updateSelectInput(session, "select_7",
                        choices = names(metab_all)[10:71])
    if (dataset_page_2 == "nano_24")
      updateSelectInput(session, "select_7",
                        choices = names(nano_24)[6:31])
    if (dataset_page_2 == "biocrates_all")
      updateSelectInput(session, "select_7",
                        choices = names(biocrates_all)[6:173])
    if (dataset_page_2 == "qpcr_all")
      updateSelectInput(session, "select_7",
                        choices = names(qpcr_all[10:19]))
    if (dataset_page_2 == "urea_all")
      updateSelectInput(session, "select_7",
                        choices = names(urea_all[8:9]))
    if (input$select_5 == "proteo_cb_sig_only")
      updateSelectInput(session, "select_7",
                        choices = names(proteo_cb_sig_only[5:25]))
    if (input$select_5 == "proteo_mctx_sig_only")
      updateSelectInput(session, "select_7",
                        choices = names(proteo_mctx_sig_only[5:21]))
    if (input$select_5 == "proteo_str_sig_only")
      updateSelectInput(session, "select_7",
                        choices = names(proteo_str_sig_only[5:43]))
    if (dataset_page_2 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_7",
                        choices = names(paul_rnaseq_pval[3:819]))
    if (dataset_page_2 == "neuron_rnaseq")
      updateSelectInput(session, "select_7",
                        choices = names(neuron_rnaseq[3:47]))
    
  }),
  
  output$Page_2_description <- renderText({"Here we are looking at the combined effects of two variables from the same dataset. By selecting two different
    variables, the correlation of the two variables together in the different traits can be visualised.
    NOTE: Axes scales differ for each plot."}),
  
  page_2_df <- reactive({
    x <- get(input$select_5)
  }),
  
  #DL Control and transgenic correlations.
  
  output$CT_comparison <- renderText({
    p2.df = page_2_df()
    paste(p2.df$tissue[1], "Fisher R-Z correlation statistic between controls and transgenics.")
  }),
  
  output$two_variable_C_and_T_comparison <- renderText({
    p2.df = page_2_df()
    paste("Fisher-value = ", round(C_and_T_compcor_list[which(names(C_and_T_compcor_list) == p2.df$tissue[1])][[1]][which(C_and_T_compcor_list[which(names(C_and_T_compcor_list) == p2.df$tissue[1])][[1]][,1] == input$select_6 & C_and_T_compcor_list[which(names(C_and_T_compcor_list) == p2.df$tissue[1])][[1]][,2] == input$select_7 | C_and_T_compcor_list[which(names(C_and_T_compcor_list) == p2.df$tissue[1])][[1]][,1] == input$select_7 & C_and_T_compcor_list[which(names(C_and_T_compcor_list) == p2.df$tissue[1])][[1]][,2] == input$select_6),][7],4))
  }),
  
  output$CT_combined_dataset_compcor_plot <- renderPlot({
    p2.df = page_2_df()
    gg.corr.vis(input$select_6, input$select_7, p2.df, "status", "Control Correlation", "Transgenic Correlation")
  }),
  
  output$pdflink2 <- downloadHandler(
    filename = function(){
      page_2_df <- reactive({
        x <- get(input$select_5)
      })
      p2.df = page_2_df()
      paste(p2.df$tissue[[1]], input$select_6, "vs", input$select_7, ".pdf")
    },
    content = function(file) {
      page_2_df <- reactive({
        x <- get(input$select_5)
      })
      p2.df = page_2_df()
      pdf(file, width = 7.5, height = 5)
      gg.corr.vis(input$select_6, input$select_7, p2.df, "status", "Control Correlation", "Transgenic Correlation")
      dev.off()
    }
  ),
  
  
  ###############################################################################################################################################################  
  ###############################################################################################################################################################
  
  ############# Differential correlation statistics ##############
  observe({
    dataset_page_5 <- input$select_13
    if (dataset_page_5 == "metab_all")
      updateSelectInput(session, "select_14",
                        choices = list("Metabolites (GC-MS) Cerebellum" = "CB_C_and_T_cors_compcor",
                                       "Metabolites (GC-MS) Hippocampus" = "hipp_C_and_T_cors_compcor",
                                       "Metabolites (GC-MS) Liver" = "liv_C_and_T_cors_compcor",
                                       "Metabolites (GC-MS) Motor Cortex" = "mctx_C_and_T_cors_compcor"))
    if (dataset_page_5 == "nano_24")
      updateSelectInput(session, "select_14",
                        choices = list("Dorsolateral" = "DL_C_and_T_cors_compcor", 
                                       "Dorsomedial" = "DM_C_and_T_cors_compcor"))
    if (dataset_page_5 == "biocrates_all")
      updateSelectInput(session, "select_14",
                        choices = list("Metabolites (LC-MS) Cerebellum" = "CB_biocrates_C_and_T_cors_compcor",
                                       "Metabolites (LC-MS) Hippocampus" = "plasma_biocrates_C_and_T_cors_compcor",
                                       "Metabolites (LC-MS) Liver" = "liv_biocrates_C_and_T_cors_compcor",
                                       "Metabolites (LC-MS) Motor Cortex" = "mctx_biocrates_C_and_T_cors_compcor"))
    #if(dataset_page_5 == "qpcr_all")
    #  updateSelectInput(session, "select_14",
    #                    choices = list("Motor Cortex QPCR" = "qpcr_mctx_C_and_T_cors_compcor",
    #                                   "Anti-Striatum QPCR" = "qpcr_antistr_C_and_T_cors_compcor",
    #                                   "Foot Plate QPCR" = "qpcr_fp_C_and_T_cors_compcor"))
    #if (dataset_page_5 == "urea_all")
    #  updateSelectInput(session, "select_14",
    #                    choices = list("Urea Serum" = "urea_serum_C_and_T_cors_compcor",
    #                                   "Urea Urine" = "urea_urine_C_and_T_cors_compcor",
    #                                   "Urea Cerebellum" = "urea_cb_C_and_T_cors_compcor",
    #                                   "Urea Hippocampus" = "urea_hipp_C_and_T_cors_compcor",
    #                                   "Urea Motor Cortex" = "urea_mctx_C_and_T_cors_compcor",
    #                                   "Urea Striatum" = "urea_stri_C_and_T_cors_compcor",
    #                                   "Urea Bladder" = "urea_blad_C_and_T_cors_compcor",
    #                                   "Urea Heart" = "urea_heart_C_and_T_cors_compcor",
    #                                   "Urea Kidney" = "urea_kid_C_and_T_cors_compcor",
    #                                   "Urea Liver" = "urea_liv_C_and_T_cors_compcor",
    #                                   "Urea Testes" = "urea_testes_C_and_T_cors_compcor"))
    if (dataset_page_5 == "Proteomics")
      updateSelectInput(session, "select_14",
                        choices = list("Cerebellum" = "proteo_cb_sig_only_C_and_T_cors_compcor",
                                       "Motor Cortex" = "proteo_mctx_sig_only_C_and_T_cors_compcor",
                                       "Striatum" = "proteo_str_sig_only_C_and_T_cors_compcor"))
    if (dataset_page_5 == "paul_rnaseq_pval")
      updateSelectInput(session, "select_14",
                        choices = list("Striatum RNASeq data" = "paul_rnaseq_pval_C_and_T_cors_compcor"))
    if (dataset_page_5 == "neuron_rnaseq")
      updateSelectInput(session, "select_14",
                        choices = list("Striatal Neuron RNASeq data" = "neuron_rnaseq_C_and_T_cors_compcor"))
    
  }),
  page_5_df <- reactive({
    t = get(input$select_14)
  }),
  
  
  output$Top_10_list <- DT::renderDataTable({
    df_3 = page_5_df()
    na.omit(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[4]] > 0 & df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[6]] > 0,][1:50,1:7])},
    rownames = F,
    options = list(
      rowCallback = JS(
        "function(row, data) {",
        "for (i = 1; i < data.length; i++) {",
        "if (data[i] > 1000 | data[i] <1){",
        "$('td:eq('+i+')', row).html(data[i].toExponential(3));",
        "}",
        "}",
        "}"),
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#00467f', 'color': '#fff'});",
        "}")
      
      )
    ),
  
  
  
  output$text_file <- downloadHandler(
    filename = function(){
      page_5_df <- reactive({
        t = get(input$select_14)
      })
      df_3 = page_5_df()
      paste(input$select_14, ".csv", sep = "")
    },
    content = function(file) {
      page_5_df <- reactive({
        t = get(input$select_14)
      })
      df_3 = page_5_df()
      rownames(df_3) <- NULL
      #write.csv(matrix(unlist(strsplit(paste(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[4]] > 0 & df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[6]] > 0,][[1]],
      #                df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[4]] > 0 & df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[6]] > 0,][[2]],
      #                sep = ", "), split = ", ")), ncol = 2, byrow = TRUE),
      #            file)
      data_to_write = na.omit(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[4]] > 0 & df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][which(df_3[order(df_3[colnames(df_3) == colnames(df_3)[7]]),][7] < 0.05),][[6]] > 0,])
      write.csv(data_to_write, file, row.names = FALSE)
      }
  ),
  
  ###############################################################################################################################################################  
  ###############################################################################################################################################################
  
  ######## Gene Set Enrichment ########
  
  #page_7_df <- reactive({
  #  x <- get(input$select_9)
  #}),
  #output$Gene_set_table <- DT::renderDataTable({
  #  p7.df = page_7_df()
  #  p7.df[,c(1, 2, 3, 7, 4, 5, 6, 13)]}, 
  #  rownames = F, 
  #  options = list(
  #    columnDefs = list(list(
  #    targets = c(0,1,2,3),
  #    render = JS(
  #      "function(data, type, row, meta) {",
  #      "return type === 'display' && data.length > 6 ?",
  #      "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
  #      "}"))
  #    ),
  #  initComplete = JS(
  #    "function(settings, json) {",
  #    "$(this.api().table().header()).css({'background-color': '#00467f', 'color': '#fff'});",
  #    "}"),
  #  rowCallback = JS(
  #    "function(row, data) {",
  #    "for (i = 1; i < data.length; i++) {",
  #    "if (data[i] > 1000 | data[i] <1){",
  #    "$('td:eq('+i+')', row).html(data[i].toExponential(3));",
  #    "}",
  #    "}",
  #    "}")
  #  )
  #  ),
  
  #output$pdflink7 <- downloadHandler(
  #  filename = function(){
  #    page_7_df <- reactive({
  #      x <- get(input$select_9)
  #    })
  #    p7.df = page_7_df()
  #    paste(input$select_9, ".csv", sep = "")
  #  },
  #  content = function(file) {
  #    page_7_df <- reactive({
  #      x <- get(input$select_9)
  #    })
  #    p7.df = page_7_df()
  #    write.csv(p7.df, file, row.names = FALSE)
  #  }
  #),
  
  
  ###############################################################################################################################################################  
  ###############################################################################################################################################################
  
  ############ Datasets Page ############
  
  output$Datasets <- renderText({"Here, you are able to download the datasets used in this application."}),
  
  output$downloadData <- downloadHandler(
    filename = function() {
      dataset_dl <- reactive({
        x <- get(input$dataset)
      })
      data = dataset_dl()
      paste(data_set_names[[which(names(data_set_names) == input$dataset)]], ".csv", sep = "")
    },
    content = function(file) {
      dataset_dl <- reactive({
        x <- get(input$dataset)
      })
      data = dataset_dl()
      write.csv(data, file, row.names = FALSE)
    }
  ),
  
  ###############################################################################################################################################################  
  ###############################################################################################################################################################
  
  ############# About ##############
  
  #par(mfrow=c(2,4)),
  #output$BRNZ <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/BRNZ.png'))
  #  list(src = filename, width = "150px", height = "150px")
  #},
  #deleteFile = FALSE),
  
  #output$CHDI <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/CHDI.png'))
  #  list(src = filename, width = "150px", height = "100px")
  #},
  #deleteFile = FALSE),
  
  #output$UOA <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/UOA.png'))
  #  list(src = filename, width = "150px", height = "75px")
  #},
  #deleteFile = FALSE),
  
  #output$SARDI <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/SARDI.png'))
  #  list(src = filename, width = "75px", height = "150px")
  #},
  #deleteFile = FALSE),
  
  #output$UOAD <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/UniofAdelaide.png'))
  #  list(src = filename, width = "150px", height = "90px")
  #},
  #deleteFile = FALSE),
  
  #output$SAHMRI <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/SAHMRI.png'))
  #  list(src = filename, width = "150px", height = "75")
  #},
  #deleteFile = FALSE),
  
  #output$CBR <- renderImage({
  #  par(mar=c(5, 5, 5, 3))
  #  filename <- normalizePath(file.path('/Users/kiwimjg/Documents/University/Summer_Schol/Emily_code/App_code/CBR.png'))
  #  list(src = filename, width = "150px", height = "75px")
  #},
  #deleteFile = FALSE),
  
  observeEvent(input$Send, {session$sendCustomMessage(type = 'testmessage',
                                                      message = 'Thank you for you email!')
    mailR::send.mail(from = "hdsheepuoa@gmail.com",
                     to = "r.snell@auckland.ac.nz",
                     subject = paste("HDSheep contact", input$First_Name, input$Last_Name),
                     body = paste(input$First_Name, " ", input$Last_Name, "\n", input$Email, "\n", input$Message, sep = ""),
                     smtp = list(host.name="mailhost.auckland.ac.nz", ssl=T))
    })
  
  
  )})

###############################################################################################################################################################  
###############################################################################################################################################################

# Run the application 

shinyApp(ui = ui, server = server)

###sendmail("mgra576@aucklanduni.ac.nz", "mgra576@aucklanduni.ac.nz", subject = "test", msg = "test", control = list(smtpServer="mailhost.auckland.ac.nz"))
