source("renv/activate.R")

# Configurar repositorios para Bioconductor y CRAN si hay conectividad
if (requireNamespace("BiocManager", quietly = TRUE)) {
  tryCatch(
    options(repos = BiocManager::repositories()),
    error = function(e) NULL
  )
}
