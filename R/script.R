#TITLE:
#"BRAF somatic mutation contributes to intrinsic epileptogenicity
#in pediatric brain tumors"

# Para acceder a los datos de secuenciación de RNA.
library("recount3")

# Guarda en una variable los proyectos de humanos disponibles.
human_projects <- available_projects()

# Selecciona los datos que concuerdan con los buscados en el proyecto SRP153743.
proj_info <- subset(
   human_projects,
   project=="SRP153743" & project_type=="data_sources")

# Crea el objeto rse para comenzar con su manipulación y análisis.
rse_gene_SRP153743 <- create_rse(proj_info)

# Cálculo de las lecturas se crean en un assay nuevo.
assay(rse_gene_SRP153743, "counts") <- compute_read_counts(rse_gene_SRP153743)

# Se obtiene la información de los metadatos.
rse_gene_SRP153743 <- expand_sra_attributes(rse_gene_SRP153743)
colData(rse_gene_SRP153743)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP153743)))
]

# Cambia los datos a numéricos y factores.
rse_gene_SRP153743$sra_attribute.age <- as.numeric(
  rse_gene_SRP153743$sra_attribute.age
)
rse_gene_SRP153743$sra_attribute.sex <- factor(
  rse_gene_SRP153743$sra_attribute.sex
)
rse_gene_SRP153743$sra_attribute.isolate <- factor(
  rse_gene_SRP153743$sra_attribute.isolate
)

# Muestra la información de los metadatos en un resumen.
summary(as.data.frame(colData(rse_gene_SRP153743)
  [
    ,
    grepl("^sra_attribute.[age|sex|isolate]", colnames(colData(rse_gene_SRP153743)))
  ]))

# Selecciona la condición de los datos.
rse_gene_SRP153743$tumor <- factor(ifelse(rse_gene_SRP153743$sra_attribute.isolate == "Tumor", "Tumor", "Normal"))
table(rse_gene_SRP153743$tumor)

# Para encontrar la proporción de lecturas asignadas a genes.
rse_gene_SRP153743$assigned_gene_prop <- rse_gene_SRP153743$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP153743$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP153743$assigned_gene_prop)

# Muestra el resumen de la proporción de lecturas asignadas a genes,
# tanto en el caso de la condición con tumor y normal.
with(colData(rse_gene_SRP153743), tapply(assigned_gene_prop, tumor, summary))

# Antes de hacer la normalización y el filtrado de datos se guarda el objeto original.
rse_gene_SRP153743_unfiltered <- rse_gene_SRP153743

# Se realiza la visualización de la proporción de lecturas asignadas a genes.
# Ajuste de margen para visualizar la imagen.
par(mar = rep(2, 4))
hist(rse_gene_SRP153743$assigned_gene_prop)

# Se filtran los datos para eliminar los genes con baja expresión.
table(rse_gene_SRP153743$assigned_gene_prop < 0.5)

# Guardar los genes con alta expresión.
rse_gene_SRP153743 <- rse_gene_SRP153743[, rse_gene_SRP153743$assigned_gene_prop >= 0.5]

# Se calcula la proporción de genes que se eliminaron.
gene_means <- rowMeans(assay(rse_gene_SRP153743, "counts"))
# Resumen del objeto anterior.
summary(gene_means)

# Se eliminan los genes con baja expresión.
rse_gene_SRP153743 <- rse_gene_SRP153743[gene_means > 0.1, ]

# Porcentaje de genes eliminados.
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

# Exploración de datos: Gráfica
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

#-------------------------------------------------------------------------------
#Verificar modelo:
# Para visualizar el modelo:
# Data frame para la visualización.
sampleData <- data.frame(
  condition = colData(rse_gene_SRP153743)$tumor,
  sex = colData(rse_gene_SRP153743)$sra_attribute.sex,
  age= colData(rse_gene_SRP153743)$sra_attribute.age,
  assigned_gene_prop = colData(rse_gene_SRP153743)$assigned_gene_prop
)
# Creación de las imágenes
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ condition + age + sex + assigned_gene_prop,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist=vd$plotlist, ncol=1)
#--------------------------------------------------------------------------------

# Usar limma para el análisis de expresión diferencial:
library("limma")
#Ajustar el tamaño de la imagen
par(mar = rep(2, 4))
# Crear el objeto DGEListSe produce la gráfica para visualizar los datos.
vGene <- voom(dge, mod, plot=TRUE)

# Utiliza el modelo bayesiano para el análisis de expresión diferencial.
# Además ajusta el modelo lineal.
eb_results <- eBayes(lmFit(vGene))

# Selecciona los genes con expresión diferencial significativa.
de_results <- topTable(
  eb_results,
  coef=2,
  number=nrow(rse_gene_SRP153743),
  sort.by="none"
)

# Pequeña exploración y visualización de los datos.
dim(de_results)
head(de_results)

# Selección de genes con una expresión diferencial significativa
table(de_results$adj.P.Val < 0.05)

# Visualización de datos
plotMA(eb_results, coef=2)

#volcano:
volcanoplot(eb_results, coef=2, highlight=4, names=de_results$gene_name)

#-------------------------------------------------------------------------------

#Visualización de genes

# Extraer los 50 genes con expresión más significativa
exprs_heatmap <- vGene$E[rank(de_results$P.Value) <= 50, ]

# Data Frame
df <- as.data.frame(colData(rse_gene_SRP153743)[, c("tumor", "sra_attribute.age", "sra_attribute.sex")])

# Nombre de las columnas del data frame
colnames(df) <- c("Condition", "AgeGroup", "Sex")

# Heatmap
#Si no se encuentra instalado el paquete, correr la siguiente línea:
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

#PCA - con la condición: Tumor / Normal
col.group <- as.factor(df$Condition)
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")

col.group <- as.character(col.group)

plotMDS(vGene$E, labels=df$Condition, col=col.group)

# PCA - sexo
col.sex <- df$Sex
levels(col.sex) <- brewer.pal(nlevels(col.sex), "Dark2")
col.sex <- as.character(col.sex)
plotMDS(vGene$E, labels=df$Sex, col=col.sex)

# PCA - edad
col.age <- df$AgeGroup
levels(col.age) <- brewer.pal(nlevels(col.age), "Set2")
col.age <- as.character(col.age)
plotMDS(vGene$E, labels=df$AgeGroup, col=col.age)
