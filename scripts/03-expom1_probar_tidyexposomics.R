# ==============================================================================
# Script: 03-expom1_probar_tidyexposomics.R
# Rama: f2-Tidyexposomics-tests
# Objetivo: Pruebas validadas contra la API exacta exportada por tidyexposomics.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyexposomics)
  library(MultiAssayExperiment)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ------------------------------------------------------------------------------
# FASE 1: CONFIGURACION Y FUNCION DE AUDITORIA
# ------------------------------------------------------------------------------
out_dir <- "results/expom_1/tidyexposomics_tests"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

test_audit <- tibble(
  Fase = character(),
  Funcion = character(),
  Estado = character(),
  Mensaje = character()
)

registrar_test <- function(fase, funcion, expr) {
  res <- tryCatch(
    {
      val <- expr
      test_audit <<- add_row(test_audit, Fase = fase, Funcion = funcion, Estado = "PASS", Mensaje = "OK")
      val
    },
    error = function(e) {
      test_audit <<- add_row(test_audit, Fase = fase, Funcion = funcion, Estado = "FAIL", Mensaje = e$message)
      NULL
    }
  )
  return(res)
}

# ------------------------------------------------------------------------------
# FASE 2: CARGA DEL OBJETO BASE (expom_1)
# ------------------------------------------------------------------------------
message("\n--- FASE 2: Cargar expom_1 ---")
load("data/raw/expom_1.RData")


# ------------------------------------------------------------------------------
# FASE 3: EXPLORACION "TIDY"
# ------------------------------------------------------------------------------
message("\n--- FASE 3: Pruebas de Manipulación Tidy ---")

registrar_test("Tidy", "pivot_sample", {
  df <- pivot_sample(expom_1)
  write_csv(head(df, 50), file.path(out_dir, "03_tidy_sample.csv"))
})

registrar_test("Tidy", "pivot_feature", {
  df <- pivot_feature(expom_1)
  write_csv(head(df, 50), file.path(out_dir, "03_tidy_feature.csv"))
})

registrar_test("Tidy", "pivot_exp", {
  df <- pivot_exp(expom_1, exp_name = "exposures")
  write_csv(head(df, 50), file.path(out_dir, "03_tidy_exp.csv"))
})

# ------------------------------------------------------------------------------
# FASE 4: VISUALIZACION DE RESULTADOS PRECOMPUTADOS (Módulos en TRUE)
# ------------------------------------------------------------------------------
message("\n--- FASE 4: Visualización ---")

# Asociación básica
registrar_test("Viz", "plot_association", {
  p <- plot_association(expom_1)
  ggsave(file.path(out_dir, "04_plot_association.png"), plot = p, width = 10, height = 13, scale = 1.1, dpi = 300)
})

# Análisis Diferencial (Volcano)
registrar_test("Viz", "plot_volcano", {
  p <- plot_volcano(expom_1)
  ggsave(file.path(out_dir, "04_plot_volcano.png"), plot = p, width = 8, height = 6)
})

# Integración DIABLO (Probando pasando el nombre de la ómica)
registrar_test("Viz", "plot_factor_summary", {
  # Si falla, el error exacto se registrará en la auditoría
  p <- plot_factor_summary(expom_1)
  ggsave(file.path(out_dir, "04_plot_factor_summary.png"), plot = p, width = 8, height = 6)
})

# Redes de co-expresión
registrar_test("Viz", "run_and_plot_network", {
  # 1. Computar la red de co-expresión/asociación para los DEGs
  expom_1 <- run_create_network(expom_1, feature_type = "degs")

  # 2. Generar el gráfico
  p_net <- plot_network(expom_1)

  # 3. Guardarlo (usamos 10x8 para que los nodos no queden apretados)
  ggsave(file.path(out_dir, "04_plot_network.png"), plot = p_net, width = 10, height = 8, dpi = 300)
})

# Enriquecimiento funcional
registrar_test("Viz", "plot_enrichment", {
  p <- plot_enrichment(expom_1)
  ggsave(file.path(out_dir, "04_plot_enrichment.png"), plot = p, width = 8, height = 6)
})

# ------------------------------------------------------------------------------
# FASE 5: EJECUCION DE NUEVAS PRUEBAS (De las funciones exportadas)
# ------------------------------------------------------------------------------
message("\n--- FASE 5: Ejecución de Nuevos Análisis ---")

# PCA del exposoma
registrar_test("Nuevos", "run_pca + plot_pca", {
  res_pca <- run_pca(expom_1)
  p <- plot_pca(res_pca)
  ggsave(file.path(out_dir, "05_plot_pca.png"), plot = p, width = 12, height = 10, scale = 1.1, dpi = 300)
})

# Normalidad
registrar_test("Nuevos", "run_normality_check", {
  res_norm <- run_normality_check(expom_1)
  p <- plot_normality_summary(res_norm)
  ggsave(file.path(out_dir, "05_plot_normality.png"), plot = p, width = 8, height = 6)
})

# Imputación
registrar_test("Nuevos", "run_impute_missing", {
  run_impute_missing(expom_1)
})

# ------------------------------------------------------------------------------
# FASE 6: AUDITORIA FINAL
# ------------------------------------------------------------------------------
message("\n--- FASE 6: Resumen y Auditoría Final ---")
audit_csv <- file.path(out_dir, "03_function_test_summary.csv")
write_csv(test_audit, audit_csv)
print(as.data.frame(test_audit))
