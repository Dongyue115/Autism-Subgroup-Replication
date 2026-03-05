rm(list = ls())

library(data.table)
library(fgsea)
library(readxl)
library(ggplot2)


e_Res <- read.csv("PLSR_result/output_gene_order_fromGLM/S1_PLS2_geneWeights.csv",
                   header = FALSE, col.names = c("Gene_symbol", "Gene_id", "Weight"))


geneList <- de_Res$Weight
names(geneList) <- de_Res$Gene_symbol


#-------------------ASD related gene sets------------------------
geneset1 <- read_excel("ASD_related genesets/geneset_004_denovo.xlsx")$gene
geneset2 <- read_excel("ASD_related genesets/geneset_003_filter1143gene.xlsx", sheet=2)$Symbol
geneset3 <- read_excel("ASD_related genesets/geneset_003_filter1143gene.xlsx", sheet=3)$Symbol
geneset4 <- read_excel("ASD_related genesets/geneset_006_105genes.xlsx")$GeneSymbol
geneset5 <- read_excel("ASD_related genesets/geneset_001.xlsx")$gene
geneset6 <- read_excel("ASD_related genesets/geneset_010_SFARI-Gene_genes_01-13-2025release_04-01-2025export.xlsx", sheet = 2)$gene_symbol
geneset7 <- read_excel("ASD_related genesets/geneset_010_SFARI-Gene_genes_01-13-2025release_04-01-2025export.xlsx", sheet = 3)$gene_symbol
geneset8 <- read_excel("ASD_related genesets/geneset_005_SPARK.xlsx")$Gene_symbol
pathway <- list(
  "ASD RDNV" = geneset1,
  "ASD Upregulated" = geneset2,
  "ASD downregulated" = geneset3,
  "FMRP-interacting" = geneset4,
  "ASD Grove et al." = geneset5,
  "Syndromic (SFARI)" = geneset6,
  "ASD SFARI" = geneset7,
  "ASD SPARK" = geneset8
)

#------------------Nonpsychiatric gene sets--------------------
# geneset1 <- read_excel("Nonpsychiatric genesets/geneset_009_psoriasis.xlsx")$gene
# geneset2 <- read_excel("Nonpsychiatric genesets/geneset_014_MSA.xlsx")$Gene_Symbol
# geneset3 <- read_excel("Nonpsychiatric genesets/geneset_012_heart_disease.xlsx")$Gene
# geneset4 <- read_excel("Nonpsychiatric genesets/geneset_014_dementia.xlsx")$Gene_Symbol
# geneset5 <- read_excel("Nonpsychiatric genesets/geneset_014_osteoarthritis.xlsx")$Gene_Symbol
# geneset6 <- read_excel("Nonpsychiatric genesets/geneset_014_CNSautoimmune.xlsx")$Gene_Symbol
# 
# 
# pathway <- list(
#   "Psoriasis" = geneset1,
#   "MSA" = geneset2,
#   "Heart disease" = geneset3,
#   "Dementia" = geneset4,
#   "Osteoarthritis" = geneset5,
#   "CNS autoimmune" = geneset6)

#------------------Other neuropsychiatric gene sets--------------------
# geneset1 <- read_excel("OtherNeuropsychiatric genesets/geneset_008_VocalLearning.xlsx")$Gene_symbol
# geneset2 <- read_excel("OtherNeuropsychiatric genesets/geneset_011_Schizophrenia.xlsx")$GeneSymbol
# geneset3 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_ADHD.xlsx")$gene
# geneset4 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_major depressive disorder.xlsx")$gene
# geneset5 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_Tourette.xlsx")$gene
# geneset6 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_generalized anxiety disorder.xlsx")$gene
# geneset7 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_antisocial personality disorder.xlsx")$gene
# geneset8 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_conduct disorder.xlsx")$gene
# geneset9 <- read_excel("OtherNeuropsychiatric genesets/geneset_007_ID.xlsx")$gene
# geneset10 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_aphasia.xlsx")$gene
# geneset11 <- read_excel("OtherNeuropsychiatric genesets/geneset_013_neurotic disorder.xlsx")$gene
# 
# pathway <- list(
#   "Vocal learning" = geneset1,
#   "Schizophrenia" = geneset2,
#   "ADHD" = geneset3,
#   "MDD" = geneset4,
#   "Tourette's" = geneset5,
#   "GAD" = geneset6,
#   "Antisocial personality disorder" = geneset7,
#   "Conduct disorder" = geneset8,
#   "ID" = geneset9,
#   "Aphasia" = geneset10,
#   "Neurotic disorder" = geneset11)


#------------------Synaptic signaling gene sets--------------------
# geneset1 <- read_excel("SynapticSignaling genesets/BP_cAMP-mediated signaling.xlsx",sheet = 1)$gene_unique
# geneset2 <- read_excel("SynapticSignaling genesets/BP_Chemical synaptic transmission.xlsx",sheet = 1)$gene_unique
# geneset3 <- read_excel("SynapticSignaling genesets/BP_G protein-coupled receptor signaling pathway.xlsx",sheet = 1)$gene_unique
# geneset4 <- read_excel("SynapticSignaling genesets/BP_Glutamate receptor signaling pathway.xlsx",sheet = 1)$gene_unique
# geneset5 <- read_excel("SynapticSignaling genesets/BP_Neurotransmitter secretion.xlsx",sheet = 1)$gene_unique
# geneset6 <- read_excel("SynapticSignaling genesets/BP_Regulation of membrane potential.xlsx",sheet = 1)$gene_unique
# geneset7 <- read_excel("SynapticSignaling genesets/BP_Synapse organization.xlsx",sheet = 1)$gene_unique
# geneset8 <- read_excel("SynapticSignaling genesets/CC_Ion channel complex.xlsx",sheet = 1)$gene_unique
# geneset9 <- read_excel("SynapticSignaling genesets/CC_Postsynaptic membrane.xlsx",sheet = 1)$gene_unique
# geneset10 <- read_excel("SynapticSignaling genesets/CC_Presynaptic active zone.xlsx",sheet = 1)$gene_unique
# geneset11 <- read_excel("SynapticSignaling genesets/CC_secretory vesicle.xlsx",sheet = 1)$gene_unique
# 
# 
# pathway <- list(
#   "cAMP-mediated signaling" = geneset1,
#   "chemical synaptic transmission" = geneset2,
#   "G protein-coupled receptor signaling pathway" = geneset3,
#   "glutamate receptor signaling pathway" = geneset4,
#   "neurotransmitter secretion" = geneset5,
#   "regulation of membrane potential" = geneset6,
#   "synapse organization" = geneset7,
#   "Ion channel complex" = geneset8,
#   "Presynaptic active zone" = geneset10,
#   "Postsynaptic membrane" = geneset9,
#   "Secretory vesicle" = geneset11)

#------------------Immune signaling gene sets--------------------
# geneset1 <- read_excel("ImmuneSignaling genesets/MF_Cytokine receptor binding.xlsx",sheet = 1)$gene_unique
# geneset2 <- read_excel("ImmuneSignaling genesets/MF_Cytokine activity.xlsx",sheet = 1)$gene_unique
# geneset3 <- read_excel("ImmuneSignaling genesets/MF_Cytokine binding.xlsx",sheet = 1)$gene_unique
# geneset4 <- read_excel("ImmuneSignaling genesets/MF_Cytokine receptor activity.xlsx",sheet = 1)$gene_unique
# geneset5 <- read_excel("ImmuneSignaling genesets/MF_Protease binding.xlsx",sheet = 1)$gene_unique
# geneset6 <- read_excel("ImmuneSignaling genesets/BP_Leukocyte chemotaxis.xlsx",sheet = 1)$gene_unique
# geneset7 <- read_excel("ImmuneSignaling genesets/BP_Response to toxic substance.xlsx",sheet = 1)$gene_unique
# geneset8 <- read_excel("ImmuneSignaling genesets/BP_Cytokine production.xlsx",sheet = 1)$gene_unique
# geneset9 <- read_excel("ImmuneSignaling genesets/BP_T cell receptor signaling pathway.xlsx",sheet = 1)$gene_unique
# geneset10 <- read_excel("ImmuneSignaling genesets/BP_Response to cytokine.xlsx",sheet = 1)$gene_unique
# geneset11 <- read_excel("ImmuneSignaling genesets/BP_Immune response transduction.xlsx",sheet = 1)$gene_unique
# geneset12 <- read_excel("ImmuneSignaling genesets/BP_Immune response.xlsx",sheet = 1)$gene_unique
# 
# pathway <- list(
#   "Cytokine receptor binding" = geneset1,
#   "Cytokine activity" = geneset2,
#   "Cytokine binding" = geneset3,
#   "Cytokine receptor activity" = geneset4,
#   "Protease binding" = geneset5,
#   "Leukocyte chemotaxis" = geneset6,
#   "Response to toxic substance" = geneset7,
#   "Cytokine production" = geneset8,
#   "T cell receptor signaling pathway" = geneset9,
#   "Response to cytokine" = geneset10,
#   "Immune response transduction" = geneset11,
#   "Immune response" = geneset12)

#------------------Protein translation gene sets--------------------
# geneset1 <- read_excel("ProteinTranslation genesets/BP_ribosomal subunit export from nucleus.xlsx",sheet = 1)$gene_unique
# geneset2 <- read_excel("ProteinTranslation genesets/CC_Ribosomal subunit.xlsx",sheet = 1)$gene_unique
# geneset3 <- read_excel("ProteinTranslation genesets/CC_Ribosome.xlsx",sheet = 1)$gene_unique
# geneset4 <- read_excel("ProteinTranslation genesets/CC_Cytosolic ribosome.xlsx",sheet = 1)$gene_unique
# geneset5 <- read_excel("ProteinTranslation genesets/BP_Cytoplasmic translation.xlsx",sheet = 1)$gene_unique
# 
# pathway <- list(
#   "Ribosomal structure" = geneset1,
#    "Ribosomal subunit" = geneset2,
#    "Ribosome" = geneset3,
#    "Cytosolic ribosome" = geneset4,
#    "Cytoplasmic translation" = geneset5)

fgseaRes <- fgsea(pathways = pathway, 
                  stats = geneList,
                  minSize = 2,
                  maxSize = 6000)


library(writexl)
write_xlsx(fgseaRes, "GSEA_result/ASD_related/s1_PLS2.xlsx")
