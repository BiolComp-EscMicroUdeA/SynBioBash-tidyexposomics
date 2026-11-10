# Guía de Alineación y Preguntas Clave para el Proyecto (scRNA-seq + Exposómica)

**Proyecto:** Disección y aplicación de `tidyexposomics` en datos de transcriptómica de célula única  
**Duración:** 2 meses (~15 horas/semana | Total: ~120 horas netas)  
**Objetivo del documento:** Definir el alcance técnico, el punto de inicio de los datos, la arquitectura computacional y los estándares de reproducibilidad para el repositorio de GitHub.

---

## 1. Escenario Upstream: Si el proyecto arranca desde FASTQ (Preprocesamiento crudo)

*En caso de que el punto de partida sean lecturas crudas, se implementará un pipeline automatizado, auditable y modular con **Snakemake** en Linux antes de pasar a R.*

1. **Tecnología de secuenciación y química de librería:**
   - ¿Qué plataforma y versión química se utilizó? (ej. 10x Genomics Chromium 3' v2/v3/v4, 5', o tecnología basada en placas como Smart-seq2)?
   - Si es 10x Genomics, ¿las lecturas están pareadas estándar (Read 1: Cell barcode + UMI; Read 2: Inserto de cDNA)?
2. **Organismo de referencia y anotación genómica:**
   - ¿Cuál es el genoma y ensamble requerido (ej. *Homo sapiens* GRCh38 / GENCODE v..., *Mus musculus* GRCm39)?
   - ¿Se debe considerar pre-mRNA / intrones (datos de *single-nucleus RNA-seq* / snRNA-seq) al momento de construir el índice de referencia?
3. **Herramienta preferida para cuantificación celular:**
   - ¿Existe preferencia dentro del grupo por alguna herramienta en particular para generar las matrices de conteos?
     - `CellRanger` (estándar de 10x, intensivo en cómputo/RAM)
     - `STARsolo` (rápido, flexible, produce salidas idénticas a CellRanger)
     - Pseudoalineadores tipo `kallisto | bustools` o `alevin-fry` (ultrarrápidos y ligeros en almacenamiento)
4. **Infraestructura y capacidad de cómputo disponible:**
   - ¿El preprocesamiento se ejecutará en un clúster local (HPC con SLURM/SGE) o en servidor dedicado? ¿Cuánta memoria RAM y CPUs están asignados para este paso?
5. **Estrategia frente a efectos de lote (Batch effects) en secuenciación:**
   - ¿Las muestras provienen de un único batch de secuenciación/célula de flujo o hubo múltiples rondas de preparación de librerías que requieran una variable de bloqueo/batch específica?

---

## 2. Nivel de Entrada de Datos y Capa de Exposoma

1. **Estado actual de los datos biológicos:**
   - Si no se arranca desde FASTQ, ¿en qué formato se transferirán las matrices de conteo? (ej. carpetas `filtered_feature_bc_matrix`, archivos `.h5`, objetos `.rds` de `SingleCellExperiment` / `Seurat`).
2. **Naturaleza y estructura del exposoma:**
   - ¿Qué tipo de exposiciones se analizarán? (químicos/metales pesados, contaminantes del aire, factores socioambientales, índices dietéticos, etc.).
   - ¿Cómo vienen estructurados los metadatos? (tabla tabular CSV/TSV con límites de detección LOD, covariables clínicas como edad, sexo, IMC, hábitos de vida).
   - ¿Las exposiciones son transversales (un punto en el tiempo por paciente) o longitudinales?

---

## 3. Integración Biológica con `tidyexposomics`

*`tidyexposomics` está diseñado sobre el ecosistema Bioconductor (`MultiAssayExperiment`, `tidySummarizedExperiment`, `tidybulk`) para vincular exposoma con capas multi-ómicas.*

1. **Resolución de la integración (Pseudobulk vs Single-Cell puro):**
   - Dado que los datos de exposición suelen medirse por donante/muestra y no por célula aislada: ¿la estrategia principal será agregar las lecturas mediante **pseudobulk por tipo celular / linaje** antes de los modelos de asociación para evitar la pseudorreplicación y el inflado de tasas de error tipo I?
   - ¿Se contempla además explorar marcadores célula a célula en subpoblaciones específicas?
2. **Modelado estadístico y contrastes:**
   - ¿Se busca un enfoque de cribado exploratorio tipo ExWAS (asociación una a una de cada exposición contra el transcriptoma)?
   - ¿O un análisis multivariado conjunto (integración de mezclas de exposición vs firmas transcriptómicas)?

---

## 4. Alcance y Entregables del Proyecto (2 Meses)

1. **Objetivo prioritario para el laboratorio:**
   - **Opción A (Metodológico/Benchmarking):** Evaluar el paquete `tidyexposomics`, sus fortalezas, limitaciones frente a otras alternativas (ej. `rexposome`/`omicRexposome`), documentando bugs o creando una viñeta/workflow de referencia.
   - **Opción B (Descubrimiento biológico):** Responder una hipótesis biológica concreta sobre la cohorte e identificar biomarcadores clave modulados por el exposoma.
2. **Productos esperados al final de las 8 semanas:**
   - Repositorio en GitHub 100% reproducible.
   - Reporte dinámico (HTML interactivo vía Quarto / R Markdown / `workflowr`).
   - Tablas de resultados (genes diferencialmente expresados, vías enriquecidas) y figuras listas para publicación.

---

## 5. Gobernanza del Repositorio de GitHub y Reproducibilidad

1. **Estándares del repositorio:**
   - ¿El laboratorio tiene directrices de ramas (`main`/`develop`), estilo de código (ej. *tidyverse style guide*) o plantillas de issues/PRs?
   - Formato preferido para reportes: ¿**Quarto**, **R Markdown**, o un sitio reproducible con **`workflowr`**?
2. **Gestión de dependencias:**
   - Se planea congelar el entorno de R utilizando `renv` (`renv.lock`) configurado con los repositorios oficiales de Bioconductor. ¿El PI o sus colaboradores usan Docker/Apptainer o prefieren restauración nativa con `renv::restore()`?
3. **Privacidad de los datos:**
   - ¿Qué nivel de confidencialidad tienen los metadatos clínicos/exposoma? (confirmar qué archivos deben estar estrictamente excluidos en `.gitignore` para no exponer información protegida en el repositorio público/privado).



# Estructura del repositorio (archivos que deben existir)

```
scRNAseq/
├── config/                                 # Parámetros y metadatos
│   ├── config.yaml                         # Parámetros del flujo (rutas, herramientas, hilos)
│   ├── samples.tsv                         # Tabla de muestras (sample_id, fastq_1, fastq_2, batch, etc.)
│   ├── units.tsv                           # Información por lane o corrida técnica (si aplica)
│   └── exposomics_metadata.tsv             # Variables ambientales, clínicas o exposómicas asociadas
│
├── workflow/                               # Flujo de trabajo Snakemake (FASTQ -> Conteos)
│   ├── Snakefile                           # Archivo principal de Snakemake
│   ├── rules/                              # Reglas modulares de Snakemake
│   │   ├── qc.smk                          # FastQC, MultiQC de lecturas crudas
│   │   ├── trimming.smk                    # Fastp / Trimmomatic (si es necesario)
│   │   ├── quantification.smk              # Alinear/cuantificar (STARsolo / CellRanger / Alevin-fry)
│   │   └── count_matrix.smk                # Filtrado inicial y exportación de matriz
│   ├── envs/                               # Definición de entornos Conda / Mamba
│   │   ├── qc_env.yaml
│   │   ├── aligner_env.yaml
│   │   └── r_env.yaml
│   └── scripts/                            # Scripts auxiliares para el workflow
│       └── make_count_matrix.R             # Genera o estandariza el formato final de matriz
│
├── analysis/                               # Análisis downstream en R (tidyexposomics)
│   ├── 00_import_and_qc.R                  # Carga matriz de Snakemake + QC inicial de células
│   ├── 01_prep_exposomics.R                # Integración con variables de exposición / metadatos
│   ├── 02_tidyexposomics_eda.R             # Análisis exploratorio (distribuciones, correlaciones)
│   ├── 03_association_testing.R            # Modelos de asociación (exp-transcriptoma / firmas)
│   └── notebooks/                          # Quarto / RMarkdown para reportes interactivos
│       └── exposomics_report.qmd
│
├── resources/                              # Referencias genómicas (o enlaces simbólicos)
│   ├── genome.fa
│   ├── annotation.gtf
│   └── index/                              # Índice del alineador (STAR, Salmon/Kallisto)
│
├── data/                                   # Datos del proyecto (usualmente no se trackean en Git)
│   ├── raw/                                # FASTQ crudos (solo lectura)
│   │   ├── sampleA_R1.fastq.gz
│   │   └── sampleA_R2.fastq.gz
│   └── exposomics/                         # Datos crudos de exposoma (mediciones, encuestas)
│
├── results/                                # Salidas generadas por Snakemake
│   ├── qc/                                 # Reportes MultiQC pre/post
│   ├── aligned/                            # BAMs/salidas intermedias
│   └── counts/                             # Matrices de conteos finales (consumidas por R)
│       ├── raw_counts/                     # Matriz cruda (ej. matrix.mtx, barcodes.tsv, features.tsv)
│       └── filtered_counts/                # Opcional: filtrado libre de células vacías (emptyDrops)
│
├── output/                                 # Resultados del análisis exposómico (análisis R)
│   ├── figures/                            # Gráficos (PCA, UMAP, volcanos, plots exposómicos)
│   ├── tables/                             # Tablas de significancia, genes DE, métricas
│   └── reports/                            # Renderizados HTML/PDF de Quarto/RMarkdown
│
├── .gitignore
├── .Rprofile
└── README.md                               # Documentación general de uso y ejecución
