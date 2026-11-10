# Datos Crudos (Raw Data)

Esta carpeta contiene los datos de entrada originales e inmutables. 

> **Regla de oro de reproducibilidad:** Los archivos en este directorio **nunca deben modificarse manualmente**. Todo procesamiento debe realizarse mediante scripts reproducibles que guarden sus salidas en `data/processed/`.

## Archivos esperados
- Datos de exposición (exposoma) en formatos tabulares (CSV, TSV).
- Matrices de expresión génica / conteos transcriptómicos (o conjuntos de datos de prueba integrados).
- Metadatos de muestras y covariables clínicas.

*Nota: Por razones de privacidad y tamaño, los datos crudos no se suben al repositorio de control de versiones (ver `.gitignore`).*

