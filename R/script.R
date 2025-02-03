#TITLE: 	BRAF somatic mutation contributes to intrinsic epileptogenicity in pediatric brain tumors

library("recount3")

human_projects <- available_projects()

proj_info <- subset(
   human_projects,
   project=="SRP153743" & project_type=="data_sources")

rse_gene_SRP153743 <- create_rse(proj_info)

assay(rse_gene_SRP153743, "counts") <- compute_read_counts(rse_gene_SRP153743)

rse_gene_SRP153743 <- expand_sra_attributes(rse_gene_SRP153743)
colData(rse_gene_SRP153743)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP153743)))
]

rse_gene_SRP153743$sra_attribute.age <- as.numeric(
  rse_gene_SRP153743$sra_attribute.age
)
rse_gene_SRP153743$sra_attribute.sex <- factor(
  rse_gene_SRP153743$sra_attribute.sex
)
rse_gene_SRP153743$sra_attribute.isolate <- factor(
  rse_gene_SRP153743$sra_attribute.isolate
)

summary(as.data.frame(colData(rse_gene_SRP153743)
  [
    ,
    grepl("^sra_attribute.[age|sex|isolate]", colnames(colData(rse_gene_SRP153743)))
  ]))

rse_gene_SRP153743$tumor <- factor(ifelse(rse_gene_SRP153743$sra_attribute.isolate == "Tumor", "Tumor", "Normal"))
table(rse_gene_SRP153743$tumor)

rse_gene_SRP153743$assigned_gene_prop <- rse_gene_SRP153743$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP153743$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP153743$assigned_gene_prop)

with(colData(rse_gene_SRP153743), tapply(assigned_gene_prop, tumor, summary))

rse_gene_SRP153743_unfiltered <- rse_gene_SRP153743

hist(rse_gene_SRP153743$assigned_gene_prop)

table(rse_gene_SRP153743$assigned_gene_prop < 0.5)

rse_gene_SRP153743 <- rse_gene_SRP153743[, rse_gene_SRP153743$assigned_gene_prop >= 0.5]

gene_means <- rowMeans(assay(rse_gene_SRP153743, "counts"))

summary(gene_means)

rse_gene_SRP153743 <- rse_gene_SRP153743[gene_means > 0.1, ]

round(nrow(rse_gene_SRP153743) / nrow(rse_gene_SRP153743_unfiltered) * 100, 2)

# Normalización de los datos

