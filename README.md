# cancer_genomics_class

Repositorio de la clase de análisis de genomas de cáncer para el Magíster en Informática Médica, realizada el 13 de octubre de 2026. Reúne la presentación de la clase y los materiales para las actividades prácticas de análisis de variantes somáticas y firmas mutacionales en cáncer colorrectal.

## Organización del repositorio

- **`class/`**: contiene la [presentación de la clase en PDF](class/Sesión12_CancerGenomes-2.pdf).
- **`practico/`**: reúne los materiales para las actividades prácticas:
  - **`scripts/`**: contiene los documentos R Markdown (`.Rmd`), sus resultados en Markdown (`.md`) y HTML, y las carpetas con las figuras generadas.
  - **`files/`**: contiene los archivos de entrada: mutaciones somáticas, datos clínicos de muestras y pacientes, y duplicación de genoma completo (WGD) junto con pureza tumoral. Adicionalmente, para revisar, se encuentran datos de alteraciones del número de copias (CNA), variantes estructurales (SV) y ancestría genética.

## Análisis prácticos

- [**Análisis de variantes somáticas**](practico/scripts/analisis_variantes_somaticas.Rmd): exploración del MAF con `maftools`, análisis de multi-hit y TMB, asociación del estado mutacional con MSI y ancestría europea, y enriquecimiento Hallmark de genes asociados a MSI.
- [**Análisis de firmas mutacionales SBS**](practico/scripts/analisis_SBS.Rmd): estimación de actividades absolutas y relativas de firmas COSMIC, visualización de prevalencia y clustering de muestras, y asociación con ancestría AMR.

## Contacto

En caso de dudas, escribir a [noemi.gomez@postgrado.uoh.cl](mailto:noemi.gomez@postgrado.uoh.cl).
