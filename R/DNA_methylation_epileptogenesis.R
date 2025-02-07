# TITULO:
# "Neuronal and glial DNA methylation and gene expression changes in early epileptogenesis"

# Abrir librería
library("recount3")

# Cargar datos de expresión génica
mouse_projects <- available_projects("mouse")
# Crear objeto de rse para manipularlo con los datos de expresión génica.
rse_gene_SRP223512 <- create_rse(
  subset(
    mouse_projects,
    project == "SRP223512" & project_type == "data_sources"
  )
)

assay(rse_gene_SRP223512, "counts") <- compute_read_counts(rse_gene_SRP223512)

# Explorar objeto de rse
rse_gene_SRP223512$sra.sample_attributes[1:5]

rse_gene_SRP223512 <- expand_sra_attributes(rse_gene_SRP223512)
colData(rse_gene_SRP223512)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP223512))
  )
]

# Cambiar las dos variables a factores:
rse_gene_SRP223512$sra_attribute.treatment <- factor(
  rse_gene_SRP223512$sra_attribute.treatment
)
rse_gene_SRP223512$sra_attribute.neun <- factor(
  rse_gene_SRP223512$sra_attribute.neun
)

# Variables de interés
summary(as.data.frame(colData(rse_gene_SRP223512)
  [
    ,
    grepl("^sra_attribute.[treatment|neun]", colnames(colData(rse_gene_SRP223512)))
  ]))

rse_gene_SRP223512$treatment <- factor(ifelse(rse_gene_SRP223512$sra_attribute.treatment == "Kainate", "epilepsy", "placebo"))
table(rse_gene_SRP223512$treatment)

rse_gene_SRP223512$assigned_gene_prop <- rse_gene_SRP223512$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP223512$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP223512$assigned_gene_prop)

with(colData(rse_gene_SRP223512), tapply(assigned_gene_prop, treatment, summary))

rse_gene_SRP223512_unfiltered <- rse_gene_SRP223512

par(mar = rep(2, 4))
hist(rse_gene_SRP223512$assigned_gene_prop)

table(rse_gene_SRP223512$assigned_gene_prop < 0.3)
rse_gene_SRP223512 <- rse_gene_SRP223512[, rse_gene_SRP223512$assigned_gene_prop > 0.3]

gene_means <- rowMeans(assay(rse_gene_SRP223512, "counts"))
summary(gene_means)

rse_gene_SRP223512 <- rse_gene_SRP223512[gene_means > 0.1, ]

dim(rse_gene_SRP223512)

round(nrow(rse_gene_SRP223512) / nrow(rse_gene_SRP223512_unfiltered) * 100, 2)

#--------------------------------------------------------------------------------

# Normalización de los datos:

library("edgeR")

dge <- DGEList(
  counts = assay(rse_gene_SRP223512, "counts"),
  genes = rowData(rse_gene_SRP223512)
  )

dge <- calcNormFactors(dge)

# ------------------------------------------------------------------------------
# Exploración de los datos:

library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP223512)), aes(y=assigned_gene_prop, x=treatment)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Treatment")

# ------------------------------------------------------------------------------
# Modelo estadístico:

mod <- model.matrix(
  ~ treatment + sra_attribute.neun + assigned_gene_prop,
  data = colData(rse_gene_SRP223512)
)
colnames(mod)

# ------------------------------------------------------------------------------
# Comprobación del modelo estadístico:
# Para visualizar el modelo:
# Data frame para la visualización.
sampleData <- data.frame(
  treatment = rse_gene_SRP223512$treatment,
  sra_attribute.neun = rse_gene_SRP223512$sra_attribute.neun,
)
# Creación de las imágenes:
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ 0 + treatment + sra_attribute.neun,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist=vd$plot, nrow = 1)

#--------------------------------------------------------------------------------

# Usar limma para el análisis de expresión diferencial:
library("limma")
# Crear el objeto DGEList
# Se produce la gráfica para visualizar los datos.
vGene <- voom(dge, mod, plot=T)

# Utiliza el modelo bayesiano para el análisis de expresión diferencial
# Además, ajusta el modelo lineal
eb_results <- eBayes(lmFit(vGene))

# Selecciona los genes con expresión diferencial significativa.
de_results <- topTable(
  eb_results,
  coef=2,
  number=nrow(rse_gene_SRP223512),
  sort.by="none"
)

# Prqueña exploración y visualización de los datos.
dim(de_results)
head(de_results)

# Selección de genes con una expresión diferencial significativa
table(de_results$adj.P.Val < 0.05)

# Visualización de datos - plotMA
plotMA(eb_results, coef=2)

# Visualización de datos - Volcano:
# REVISAR!!!!!
volcanoplot(eb_results, coef=2, highlight=4, names=de_results$gene_name)

#--------------------------------------------------------------------------------

# Visualización de genes

# Extraer los 50 genes con expresión más significativa
exprs_heatmap <- vGene$E[rank(de_results$P.Value) <= 50, ]

# Data frame
df <- as.data.frame(colData(rse_gene_SRP223512)[, c("treatment", "sra_attribute.neun")])

# Nombre de las columnas del data frame
colnames(df) <- c("treatment", "Neun/marcador")

# Heatmap
# Si no se encuentra instalado el paquete, correr la siguiente línea:
# install.packages("pheatmap")
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = F,
  show_colnames = F,
  annotation_col = df,
)

library("RColorBrewer")

# PCA - con el tratamiento: Kainate/Sham
col.group <- as.factor(df$treatment)
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")

col.group <- as.character(col.group)

plotMDS(vGene$E, labels = df$treatment, col=col.group)

# PCA - con el Neun: positivo/negativo
col.neun <- df$`Neun/marcador`
levels(col.neun) <- brewer.pal(nlevels(col.neun), "Dark2")
col.neun <- as.character(col.neun)
plotMDS(vGene$E, labels = df$`Neun/marcador`, col=col.neun)
