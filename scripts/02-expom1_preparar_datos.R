# ============================================================
# Script: 02-expom1_preparar_datos.R
# Rama: p01/inspection
# Objetivo: Auditar el objeto expom_1 (MultiAssayExperiment),
#           inspeccionar resultados preexistentes embebidos,
#           analizar missingness y perfilar expresion genica.
# ============================================================


# 1. CONFIGURACION --------------------------------------------

library(tidyexposomics)
library(MultiAssayExperiment)
library(SummarizedExperiment)

data_dir <- "data/raw"
results_dir <- "results/expom_1"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat(
  "tidyexposomics:",
  as.character(packageVersion("tidyexposomics")),
  "\n"
)


# 2. CARGA Y VALIDACION DEL MAE -------------------------------

env <- new.env()

load(
  file.path(data_dir, "expom_1.RData"),
  envir = env
)

expom <- env$expom_1
rm(env)

stopifnot(inherits(expom, "MultiAssayExperiment"))

cat("\nObjeto cargado correctamente.\n")
cat("Ensayos:", names(experiments(expom)), "\n")


# 3. AUDITORIA DE RESULTADOS PREEXISTENTES ---------------------

# 3.1 Inventario de resultados embebidos

cat("\n--- EXISTING RESULTS ---\n")

result_types <- c(
  "codebook",
  "quality_control",
  "correlation",
  "association",
  "mixture_analysis",
  "differential_analysis",
  "multiomics_integration",
  "network",
  "enrichment"
)

result_inventory <- data.frame(
  result_type = result_types,
  available = vapply(
    result_types,
    function(x) {
      !is.null(extract_results(
        exposomicset = expom,
        result = x
      ))
    },
    logical(1)
  )
)

print(result_inventory)

write.csv(
  result_inventory,
  file.path(results_dir, "02_existing_results.csv"),
  row.names = FALSE
)


# 3.2 Resumen del pipeline (si existe historial)

cat("\n--- PIPELINE SUMMARY ---\n")

expom_metadata <- metadata(expom)

has_pipeline <- (
  !is.null(expom_metadata$summary) &&
  !is.null(expom_metadata$summary$steps)
)

if (has_pipeline) {

  run_pipeline_summary(
    exposomicset = expom,
    console_print = TRUE,
    diagram_print = FALSE,
    include_notes = TRUE
  )

} else {

  cat("No se encontro historial de pipeline.\n")

}


# 4. INVENTARIO DE COLDATA ------------------------------------

cat("\n--- COLDATA INVENTORY ---\n")

cd <- as.data.frame(colData(expom))

coldata_inventory <- data.frame(
  variable = names(cd),
  class = vapply(
    cd,
    function(x) class(x)[1],
    character(1)
  ),
  missing = colSums(is.na(cd)),
  unique_values = vapply(
    cd,
    function(x) length(unique(x[!is.na(x)])),
    integer(1)
  )
)

write.csv(
  coldata_inventory,
  file.path(
    results_dir,
    "02_coldata_inventory.csv"
  ),
  row.names = FALSE
)

cat("Variables en colData:", ncol(cd), "\n")
cat("Inventario exportado.\n")


# 5. CONTROL DE CALIDAD: MISSINGNESS --------------------------

cat("\n--- MISSINGNESS ANALYSIS ---\n")

# 5.1 Resumen general de missingness

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

# 5.2 Missingness por exposicion (lollipop)

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

# 5.3 Tabla de variables con valores faltantes

missing_summary <- data.frame(
  variable = names(cd),
  missing_values = colSums(is.na(cd)),
  missing_pct = colMeans(is.na(cd)) * 100
)

missing_summary <- missing_summary[
  missing_summary$missing_values > 0,
]

missing_summary <- missing_summary[
  order(-missing_summary$missing_pct),
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
  "\nVariables con valores faltantes:",
  nrow(missing_summary), "\n"
)

print(head(missing_summary, 15))

cat("\nQuality control completado.\n")


# 6. PERFILADO DE EXPRESION GENICA -----------------------------

cat("\n--- EXPRESSION SUMMARY ---\n")

gene_expression <- experiments(expom)[["Gene Expression"]]

cat("Clase:", class(gene_expression), "\n")
cat("Dimensiones:", dim(gene_expression), "\n")

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
  file.path(
    results_dir,
    "02_gene_expression_variance.csv"
  ),
  row.names = FALSE
)

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


# 7. INSPECCION DE RESULTADOS EMBEBIDOS ------------------------

cat("\n--- INSPECTING EMBEDDED RESULTS ---\n")

inspect_result <- function(result_name) {

  result <- tryCatch(
    extract_results(
      exposomicset = expom,
      result = result_name
    ),
    error = function(e) {
      cat("Error extrayendo", result_name, ":", e$message, "\n")
      return(NULL)
    }
  )

  cat("\n---", result_name, "---\n")
  cat(
    "Clase:",
    paste(class(result), collapse = ", "),
    "\n"
  )

  if (is.null(result)) {
    cat("Resultado no disponible.\n")
    return(invisible(NULL))
  }

  if (is.list(result)) {

    cat("Elementos:", length(result), "\n")
    print(head(names(result), 15))

    n_show <- min(length(result), 10)

    for (i in seq_len(n_show)) {

      x <- result[[i]]
      label <- if (!is.null(names(result))) {
        names(result)[i]
      } else {
        as.character(i)
      }

      dims_str <- if (!is.null(dim(x))) {
        paste(dim(x), collapse = " x ")
      } else {
        ""
      }

      cat(
        label, ": ",
        paste(class(x), collapse = ", "),
        " ", dims_str, "\n"
      )
    }

  } else if (is.data.frame(result)) {

    cat("Dimensiones:", dim(result), "\n")
    print(head(result, 5))

  } else {

    print(utils::str(result, max.level = 1))

  }

  invisible(result)
}

# Inspeccionar resultados disponibles
available_results <- result_inventory$result_type[
  result_inventory$available
]

for (res_name in available_results) {
  inspect_result(res_name)
}


# 8. INSPECCION DIRIGIDA DE RESULTADOS ------------------------

inspect_tree <- function(x, path = "root", depth = 0, max_depth = 3) {

  cat("\n", path, "\n", sep = "")
  cat("Clase:", paste(class(x), collapse = ", "), "\n")

  if (is.data.frame(x)) {

    cat("Dimensiones:", paste(dim(x), collapse = " x "), "\n")
    cat("Columnas:", paste(names(x), collapse = ", "), "\n")

  } else if (is.list(x) && depth < max_depth) {

    cat("Elementos:", length(x), "\n")

    for (i in seq_len(min(length(x), 8))) {

      element_name <- names(x)[i]

      if (is.null(element_name) || element_name == "") {
        element_name <- paste0("element_", i)
      }

      inspect_tree(
        x[[i]],
        path = paste(path, element_name, sep = "$"),
        depth = depth + 1,
        max_depth = max_depth
      )
    }
  }
}

# Association
association_results <- extract_results(
  exposomicset = expom,
  result = "association"
)

cat("\n========== ASSOCIATION ==========\n")
inspect_tree(association_results)

# Multiomics integration
diablo_results <- extract_results(
  exposomicset = expom,
  result = "multiomics_integration"
)

cat("\n========== DIABLO ==========\n")
inspect_tree(diablo_results)


# 8.1 DIABLO: RESUMEN DE INTEGRACION --------------------------

diablo <- diablo_results$integration_results$result
top_features <- diablo_results$top_factor_features
common_features <- diablo_results$common_top_factor_features

cat("\n--- METHOD ---\n")
print(diablo_results$integration_results$method)

cat("\n--- DIABLO MODEL ---\n")
cat(
  "Class:",
  paste(class(diablo), collapse = ", "), "\n"
)
cat("Components:", diablo$ncomp, "\n")
cat(
  "Blocks:",
  paste(names(diablo$X), collapse = ", "), "\n"
)

cat("\n--- OUTCOME CLASSES ---\n")
print(table(diablo$Y))

cat("\n--- SELECTED FEATURES ---\n")
print(table(top_features$exp_name))
print(table(top_features$component))

cat("\n--- COMMON FEATURES ---\n")
cat("Total:", nrow(common_features), "\n")
cat("Features by omics:\n")
print(table(common_features$exp_name))

cat("\n--- EXAMPLE ANNOTATED FEATURES ---\n")
print(
  head(
    common_features[, c(
      "feature", "exp_name", "component",
      "loading", "factor", "rank",
      "feature_clean"
    )],
    10
  )
)

# 9. REPRODUCIBILIDAD Y SESION --------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "02_session_info.txt")
)

cat("\n========================================\n")
cat("SCRIPT 02 COMPLETADO\n")
cat("========================================\n")
cat("Archivos generados en:", results_dir, "\n")
cat(
  " -",
  paste(
    list.files(results_dir),
    collapse = "\n - "
  ),
  "\n"
)
