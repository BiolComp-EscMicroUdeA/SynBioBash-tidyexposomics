# ==============================================================================
# 02_preparar_datos.R
# Propósito: Limpieza, formateo y armado de objetos estructurados
# (p. ej. MultiAssayExperiment / SummarizedExperiment) para tidyexposomics.
# ==============================================================================


# 1. CONFIGURACION --------------------------------------------

library(tidyexposomics)
library(MultiAssayExperiment)
library(SummarizedExperiment)

data_dir <- "data/raw"
results_dir <- "results/expom_1"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

cat(
  "tidyexposomics:",
  as.character(packageVersion("tidyexposomics")), "\n"
)


# 2. CARGAR MULTIASSAYEXPERIMENT -------------------------------

env <- new.env()

load(
  file.path(data_dir, "expom_1.RData"),
  envir = env
)

expom <- env$expom_1

stopifnot(inherits(expom, "MultiAssayExperiment"))

cat("\nObjeto cargado correctamente.\n")
cat("Ensayos:", names(experiments(expom)), "\n")


# 3. SELECCIONAR GENE EXPRESSION -------------------------------

gene_expression <- experiments(expom)[["Gene Expression"]]

cat("\n--- GENE EXPRESSION ---\n")
cat("Clase:", class(gene_expression), "\n")
cat("Dimensiones:", dim(gene_expression), "\n")


# 4. QUALITY CONTROL: MISSINGNESS ------------------------------

cat("\n--- MISSINGNESS ANALYSIS ---\n")

missing_plot <- plot_missing(
  exposomicset = expom,
  plot_type = "summary",
  threshold = 0
)

print(missing_plot)

ggplot2::ggsave(
  filename = file.path(
    results_dir,
    "02_missingness_summary.png"
  ),
  plot = missing_plot,
  width = 10,
  height = 6,
  dpi = 300
)

cat("Grafico de missingness exportado.\n")


# 4.1 IDENTIFICAR EXPOSICIONES CON MISSINGNESS --------

exposure_missing_plot <- plot_missing(
  exposomicset = expom,
  plot_type = "lollipop",
  threshold = 0,
  layers = "Exposure"
)

print(exposure_missing_plot)

ggplot2::ggsave(
  filename = file.path(
    results_dir,
    "02_exposure_missingness.png"
  ),
  plot = exposure_missing_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Conteo de NA por variable en colData
cd <- as.data.frame(colData(expom))

missing_summary <- data.frame(
  variable = names(cd),
  missing_values = colSums(is.na(cd)),
  missing_percentage = colMeans(is.na(cd)) * 100
)

missing_summary <- missing_summary[
  missing_summary$missing_values > 0,
]

missing_summary <- missing_summary[
  order(-missing_summary$missing_percentage),
]

write.csv(
  missing_summary,
  file.path(
    results_dir,
    "02_missing_values_by_variable.csv"
  ),
  row.names = FALSE
)

cat(
  "\nVariables con valores faltantes en colData:",
  nrow(missing_summary), "\n"
)

print(head(missing_summary, 15))


# 5. RESUMEN DE EXPRESION GENICA --------------------------------

expression_matrix <- assay(gene_expression)

feature_variance <- apply(
  expression_matrix,
  1,
  var,
  na.rm = TRUE
)

expression_summary <- data.frame(
  feature_id = rownames(expression_matrix),
  variance = feature_variance
)

write.csv(
  expression_summary,
  file.path(results_dir, "02_gene_expression_variance.csv"),
  row.names = FALSE
)

cat("\n--- EXPRESSION SUMMARY ---\n")
cat("Features:", nrow(expression_matrix), "\n")
cat("Muestras:", ncol(expression_matrix), "\n")
cat(
  "Features con varianza cero:",
  sum(feature_variance == 0, na.rm = TRUE), "\n"
)
cat(
  "Features con varianza NA:",
  sum(is.na(feature_variance)), "\n"
)


# 6. GUARDAR RESUMEN --------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "02_session_info.txt")
)

cat("\nExploracion finalizada.\n")
