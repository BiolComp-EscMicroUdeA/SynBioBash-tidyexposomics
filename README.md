# tidyexposomics-pruebas

Proyecto para la disección, evaluación y benchmarking del paquete **tidyexposomics** y su integración con datos transcriptómicos / exposómicos.

## Estructura del Proyecto

```text
tidyexposomics-pruebas/
├── README.md                       # Documentación general del proyecto y flujo de trabajo
├── .gitignore                      # Archivos y carpetas excluidos del control de versiones
├── tidyexposomics-pruebas.Rproj    # Archivo de configuración de proyecto RStudio
├── renv.lock                       # Bloqueo de dependencias para reproducibilidad
├── renv/                           # Configuración interna de renv
├── data/
│   ├── raw/                        # Datos originales e inmutables
│   │   └── README.md               # Descripción de origen, versión y formato de los datos crudos
│   └── processed/                  # Datos limpios y objetos serializados (.rds)
├── scripts/
│   ├── README.md                   # Descripción de la ejecución secuencial de scripts
│   ├── 01_inspeccionar_datos.R     # Exploración inicial y estructura de datos
│   ├── 02_preparar_datos.R         # Curación y construcción de contenedores (MultiAssayExperiment)
│   ├── 03_probar_tidyexposomics.R  # Pruebas de las funciones clave de tidyexposomics
│   └── 04_validar_resultados.R     # Validación y comparación de salidas
├── results/
│   ├── figures/                    # Figuras y gráficos generados
│   ├── tables/                     # Tablas de resultados (CSV/TSV)
│   └── reports/                    # Reportes compilados en HTML/PDF
└── tests/
    └── testthat/
        └── test-resultados.R       # Pruebas automatizadas de consistencia de salidas
```

## Reproducibilidad con `renv`

Para restaurar el entorno computacional exacto en otra máquina:

```r
renv::restore()
```

