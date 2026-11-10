# Scripts de Análisis y Procesamiento

Flujo de ejecución secuencial para la disección y prueba de `tidyexposomics`:

1. **`01_inspeccionar_datos.R`**: Carga y examen inicial de dimensiones, tipos de variables, datos faltantes y distribución de las matrices.
2. **`02_preparar_datos.R`**: Limpieza, alineación de identificadores de muestra y ensamblado del contenedor estructurado (ej. `MultiAssayExperiment`).
3. **`03_probar_tidyexposomics.R`**: Ejecución de las funciones principales de `tidyexposomics` (asociaciones, integración, métodos tidy).
4. **`04_validar_resultados.R`**: Verificación de robustez, chequeos cruzados de consistencia y exportación de tablas/gráficos clave.

