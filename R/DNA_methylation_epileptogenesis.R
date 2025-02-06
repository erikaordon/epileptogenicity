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
