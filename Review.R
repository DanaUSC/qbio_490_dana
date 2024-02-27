# Final R Project
# QBIO 490 Spring 2024
# Dana Souter

# Libraries
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(maftools)
library(ggplot2)
library(survival)
library(survminer)
library(data.table)

knitr::opts_knit$set(root.dir = normalizePath("/Users/danasouter/Desktop/qbio490/qbio_490_dana/analysis_data"))

# Load and prepare  data
clinic <- fread("/Users/danasouter/Desktop/qbio490/qbio_490_dana/analysis_data", data.table = FALSE)
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"


# metastatic status 
#compute survival time and death event
clinic_mod <- clinic
clinic_mod$death_event <- ifelse(clinic_mod$vital_status == "Alive", 0, 1)
clinic_mod$survival_time <- ifelse(is.na(clinic_mod$days_to_death), clinic_mod$days_to_last_follow_up, clinic_mod$days_to_death)

# Survival Analysis for Metastatic vs. Non-Metastatic
surv_object <- Surv(time = clinic_mod$survival_time, event = clinic_mod$death_event)
metastasis_fit <- survfit(surv_object ~ clinic_mod$metastatic_status, data = clinic_mod)
metastasis_plot <- ggsurvplot(metastasis_fit, pval = TRUE, ggtheme = theme_minimal(), legend = "right")

#print
metastasis_plot


# Boxplot for Survival Time by Metastatic Status
boxplot(survival_time ~ metastatic_status, data = clinic_mod, xlab = "Metastatic Status", ylab = "Survival Time")

# Load mutation data
maf_dataframe <- fread("/Users/danasouter/Desktop/qbio490/qbio_490_dana/analysis_data", data.table = FALSE)
maf_object <- read.maf(maf = maf_dataframe, clinicalData = clinic_mod, isTCGA = TRUE)

# Co-oncoplot for Metastatic vs. Non-Metastatic Patients
metastatic_ids <- clinic_mod$Tumor_Sample_Barcode[clinic_mod$metastatic_status == "Metastatic"]
non_metastatic_ids <- clinic_mod$Tumor_Sample_Barcode[clinic_mod$metastatic_status == "Non-metastatic"]
metastatic_maf <- subsetMaf(maf = maf_object, tsb = metastatic_ids)
non_metastatic_maf <- subsetMaf(maf = maf_object, tsb = non_metastatic_ids)

#print
coOncoplot(m1 = metastatic_maf, m2 = non_metastatic_maf, m1Name = "Metastatic Patients", m2Name = "Non-metastatic Patients")


# Differential Expression Analysis
dds <- DESeqDataSetFromMatrix(countData = assays(sum_exp)$counts, colData = colData(sum_exp), design = ~ metastatic_status + treatment + race + gender + vital_status)
dds <- DESeq(dds)
res <- results(dds, contrast = c("metastatic_status", "Metastatic", "Non-metastatic"))


# Volcano Plot
volcano_plot <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05 & abs(log2FoldChange) > log2(fc_threshold), "Significant", "Not significant"))) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  scale_color_manual(values = c("red", "black")) +
  labs(x = "Log2 Fold Change (Metastatic/Non-Metastatic)", y = "-Log10 Adjusted P-value")

#print
volcano_plot


