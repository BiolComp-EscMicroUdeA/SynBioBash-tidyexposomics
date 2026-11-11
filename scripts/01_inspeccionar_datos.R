# ==============================================================================
# Script: 01_inspeccionar_datos.R
# Rama: f0-Data-inspection
# Objetivo: Cargar .RData de Zenodo, auditar clases,
#           dimensiones y missing values.
# ==============================================================================


# ------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------

# Project directories
data_dir <- file.path("data", "raw")
tables_results_dir <- "results/tables"

# Expected datasets
expected_files <- c(
  "annotated_cb.RData",
  "expom_1.RData",
  "fdata.RData",
  "hca_lung.RData",
  "meta.RData",
  "omics_list.RData"
)

# Create results directory if it does not exist
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Check available files
available_files <- list.files(
  data_dir,
  pattern = "\\.RData$",
  full.names = FALSE
)

missing_files <- setdiff(expected_files, available_files)

if (length(missing_files) > 0) {
  warning(
    "Missing expected files: ",
    paste(missing_files, collapse = ", ")
  )
}

files_to_inspect <- intersect(expected_files, available_files)

if (length(files_to_inspect) == 0) {
  stop("No expected RData files were found in data/raw/.")
}

cat("Files to inspect:", length(files_to_inspect), "\n")


# ------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------

# Safely extract object dimensions
get_dimensions <- function(x) {
  dimensions <- tryCatch(
    dim(x),
    error = function(e) NULL
  )

  if (is.null(dimensions)) {
    return(NA_character_)
  }

  paste(dimensions, collapse = " x ")
}


# Generate a short preview of names
preview_names <- function(x, max_items = 10) {
  object_names <- tryCatch(
    names(x),
    error = function(e) NULL
  )

  if (is.null(object_names)) {
    return(NA_character_)
  }

  if (length(object_names) == 0) {
    return(NA_character_)
  }

  preview <- head(object_names, max_items)

  result <- paste(preview, collapse = ", ")

  if (length(object_names) > max_items) {
    result <- paste0(result, ", ...")
  }

  result
}


# ------------------------------------------------------------
# 3. INSPECTION
# ------------------------------------------------------------

inventory <- list()
structure_report <- c(
  "TIDYEXPOSOMICS DATA INSPECTION",
  "================================",
  ""
)

for (file_name in files_to_inspect) {
  file_path <- file.path(data_dir, file_name)

  cat("\nInspecting:", file_name, "\n")

  # Create an isolated environment for each file
  object_env <- new.env()

  # Load objects without polluting the global environment
  loaded_objects <- load(
    file_path,
    envir = object_env
  )

  structure_report <- c(
    structure_report,
    "",
    "----------------------------------------",
    paste("FILE:", file_name),
    "----------------------------------------",
    paste("Objects found:", length(loaded_objects)),
    ""
  )

  # Inspect every object within the file
  for (object_name in loaded_objects) {
    object <- get(
      object_name,
      envir = object_env,
      inherits = FALSE
    )

    object_class <- paste(class(object), collapse = ", ")

    object_type <- typeof(object)

    object_dimensions <- get_dimensions(object)

    object_length <- length(object)

    object_names <- preview_names(object)

    # Extract column names when applicable
    column_names <- if (is.data.frame(object)) {
      preview_names(object)
    } else {
      NA_character_
    }

    # Store inventory information
    inventory[[length(inventory) + 1]] <- data.frame(
      file = file_name,
      object = object_name,
      class = object_class,
      type = object_type,
      dimensions = object_dimensions,
      length = object_length,
      names_preview = object_names,
      columns_preview = column_names,
      stringsAsFactors = FALSE
    )

    # Print a summary to the console
    cat("  Object:", object_name, "\n")
    cat("  Class:", object_class, "\n")
    cat("  Type:", object_type, "\n")
    cat("  Dimensions:", object_dimensions, "\n")

    # Capture a compact structural representation
    object_structure <- utils::capture.output(
      utils::str(
        object,
        max.level = 2,
        list.len = 10,
        vec.len = 6,
        give.attr = FALSE
      )
    )

    # Add information to the report
    structure_report <- c(
      structure_report,
      paste("OBJECT:", object_name),
      paste("Class:", object_class),
      paste("Type:", object_type),
      paste("Dimensions:", object_dimensions),
      paste("Length:", object_length),
      "Structure:",
      object_structure,
      ""
    )
  }

  # Release memory before loading the next file
  rm(object_env)
  invisible(gc())
}


# ------------------------------------------------------------
# 4. EXPORT RESULTS
# ------------------------------------------------------------

# Combine individual records into one data frame
data_inventory <- do.call(
  rbind,
  inventory
)

# Export inventory
utils::write.csv(
  data_inventory,
  file = file.path(tables_results_dir, "01_data_inventory.csv"),
  row.names = FALSE,
  na = ""
)

# Export structural report
writeLines(
  structure_report,
  con = file.path(tables_results_dir, "01_structure_report.txt"),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 5. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n========================================\n")
cat("INSPECTION COMPLETED\n")
cat("========================================\n")

cat("Files inspected:", length(files_to_inspect), "\n")
cat("Objects identified:", nrow(data_inventory), "\n")

cat("\nInventory:\n")
print(data_inventory)

cat("\nGenerated files:\n")
cat("- results/01_data_inventory.csv\n")
cat("- results/01_structure_report.txt\n")
