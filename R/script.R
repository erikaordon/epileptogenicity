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

par(mar = rep(2, 4))
hist(rse_gene_SRP153743$assigned_gene_prop)

table(rse_gene_SRP153743$assigned_gene_prop < 0.5)

rse_gene_SRP153743 <- rse_gene_SRP153743[, rse_gene_SRP153743$assigned_gene_prop >= 0.5]

gene_means <- rowMeans(assay(rse_gene_SRP153743, "counts"))

summary(gene_means)

rse_gene_SRP153743 <- rse_gene_SRP153743[gene_means > 0.1, ]

round(nrow(rse_gene_SRP153743) / nrow(rse_gene_SRP153743_unfiltered) * 100, 2)

# ------------------------------------------------------------------------------

# Normalización de los datos

# En caso de no contar con la librería descargada:
# BiocManager::install("edgeR", update=F)
library("edgeR")

# Creación de la matriz de conteos.
dge <- DGEList(
  counts=assay(rse_gene_SRP153743, "counts"),
  genes=rowData(rse_gene_SRP153743)
)

# Normalización de los datos utilizando TMM para edgeR.
dge <- calcNormFactors(dge)

# ------------------------------------------------------------------------------

# Expresión diferencial

# Exploración de datos:
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP153743)), aes(y=assigned_gene_prop, x=tumor)) +
  geom_boxplot() +
  theme_bw(base_size=20) +
  ylab("Assigned Gene Prop") +
  xlab("Tumor Group")

# Modelo estadístico: (CHECAR)
mod <- model.matrix(~ tumor + sra_attribute.age + sra_attribute.sex + assigned_gene_prop,
  data = colData(rse_gene_SRP153743)
  )
colnames(mod)

# Usar limma para el análisis de expresión diferencial:
library("limma")
#Ajustar el tamaño de la imagen
par(mar = rep(2, 4))
# Crear el objeto DGEListSe produce la gráfica para visualizar los datos.
vGene <- voom(dge, mod, plot=TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef=2,
  number=nrow(rse_gene_SRP153743),
  sort.by="none"
)

dim(de_results)
head(de_results)

# Selección de genes con una expresión diferencial significativa
table(de_results$adj.P.Val < 0.05)

# Visualización de datos
plotMA(eb_results, coef=2)

#volcano:
volcanoplot(eb_results, coef=2, highlight=4, names=de_results$gene_name)

#Visualización de genes
