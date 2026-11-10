# ==============================================================================
# Script: 01_inspeccionar_datos.R
# Rama: f0-Data-inspection
# Objetivo: Cargar .RData de Zenodo, auditar clases,
#           dimensiones y missing values.
# ==============================================================================

library(tidyexposomics)
library(dplyr)

# 1. Definir rutas
raw_dir  <- "data/raw"
rdata_path <- list.files(raw_dir, pattern = "\\.RData$", full.names = TRUE)[1]

if (is.na(rdata_path)) {
  stop(
    "No se encontró ningún archivo .RData en 'data/raw/'.\n",
    "Ubica el archivo descargado de Zenodo allí."
  )
}

message(paste("Cargando archivo:", rdata_path))

# 2. Cargar .RData en un entorno aislado para no contaminar el global
env_data <- new.env()
loaded_objects <- load(rdata_path, envir = env_data)

message("Objetos encontrados en el .RData:")
print(loaded_objects)

#######
# 1. Extraer el objeto del entorno aislado
cb_obj <- get("annotated_cb", envir = env_data)

# 2. Verificar su clase y dimensiones
class(cb_obj)
dim(cb_obj)

# 3. Inspeccionar nombres de columnas y estructura de variables
colnames(cb_obj)
head(cb_obj)

################




# 3. Extraer metadatos de cada objeto cargado
inspeccion_df <- purrr::map_df(loaded_objects, function(obj_name) {
  obj <- get(obj_name, envir = env_data)

  # Determinar dimensiones
  dims <- if (is.data.frame(obj) || is.matrix(obj)) {
    paste(dim(obj), collapse = " x ")
  } else if (
    inherits(obj, "ExposomeSet") || inherits(obj, "MultiAssayExperiment")
  ) {
    paste(dim(obj), collapse = " x ")
  } else {
    paste("Longitud:", length(obj))
  }

  # Calcular porcentaje de NA si es data frame o matriz
  pct_na <- if (is.data.frame(obj) || is.matrix(obj)) {
    round(mean(is.na(obj)) * 100, 2)
  } else {
    NA
  }

  tibble(
    Objeto = obj_name,
    Clase = paste(class(obj), collapse = ", "),
    Dimensiones = dims,
    Porcentaje_NA = pct_na
  )
})

# 4. Mostrar resumen en consola
print(inspeccion_df)

# 5. Guardar tabla de inspección
out_table <- "results/tables/resumen_inspeccion.csv"
write_csv(inspeccion_df, out_table)
message(paste("Resumen guardado exitosamente en:", out_table))

# 6. Inspección detallada si existen tablas independientes de exposoma
if ("codebook" %in% loaded_objects) {
  cb <- get("codebook", envir = env_data)
  message("\n--- Vista previa del Codebook ---")
  print(head(cb))
}

if ("exposures" %in% loaded_objects) {
  exp_df <- get("exposures", envir = env_data)
  message("\n--- Vista previa de Exposiciones (Primeras 5 columnas) ---")
  print(head(exp_df[, 1:min(5, ncol(exp_df))]))
}