
# High-throughput transcriptome sequencing and comparative analysis of *Escherichia coli* and *Schizosaccharomyces pombe* in respiratory and fermentative growth

**Fecha:** 29/09/2025  
**Autores:** Miryam Zamora, Alondra Márquez  

---

## Descripción del proyecto

El presente proyecto tiene como propósito construir un **pipeline reproducible de RNA-seq** para ambos organismos y comparar firmas de expresión, así como analizar, mediante herramientas bioinformáticas, los resultados reportados en el estudio:

> Vichi J., Salazar E., Jiménez-Jacinto V., Olvera Rodríguez L., Grande R., Dantán-González E., Morett E., Hernández-Mendoza A. (2021).  
> *High-throughput transcriptome sequencing and comparative analysis of* Escherichia coli *and* Schizosaccharomyces pombe *in respiratory and fermentative growth*. PLOS ONE.  
> [https://doi.org/10.1371/journal.pone.0248513](https://doi.org/10.1371/journal.pone.0248513)

En dicho trabajo, se llevó a cabo la secuenciación y análisis comparativo del transcriptoma de *Escherichia coli* (procariota) y *Schizosaccharomyces pombe* (eucariota unicelular) bajo condiciones de crecimiento respiratorias y fermentativas. A partir de estos datos, se identificaron genes diferencialmente expresados y se exploraron las similitudes y diferencias en los procesos metabólicos y regulatorios entre ambos organismos.

En este proyecto se emplearán directamente los datos proporcionados por los autores para llevar a cabo el pipeline de análisis bioinformático.

---

## Objetivo general

Analizar los datos transcriptómicos de *E. coli* y *S. pombe* bajo condiciones de respiración y fermentación con el fin de:

- Identificar genes diferencialmente expresados.
- Explorar su anotación funcional.
- Realizar una comparación entre organismos que permita resaltar procesos conservados y divergentes.

---

## Objetivos específicos

1. Descargar las 12 bibliotecas FASTQ (2 organismos × 2 condiciones × 3 réplicas).
2. Descargar y versionar genomas de referencia:
   - *Schizosaccharomyces pombe* (fission yeast) ASM294v3 — `GCF_000002945.2` (reference genome).
   - *Escherichia coli* str. K-12 substr. MG1655 — `GCF_000005845.2` (reference genome).
   - Incluir secuencias en formato FASTA y anotaciones en GTF/GFF.  
     > Se descargarán automáticamente a `../../data/`.
3. Alinear por réplica y generar archivos BAM ordenados e indexados (`.bam` + `.bai`).
4. Obtener una matriz de conteos por gen.
5. Realizar análisis de expresión diferencial (DESeq2) por organismo.
6. Anotar los DEGs y comparar patrones compartidos/divergentes entre *E. coli* y *S. pombe*.

---

## Planteamiento del problema

### Contexto

La respiración y la fermentación inducen reprogramaciones transcriptómicas profundas. Comparar estas respuestas en dos organismos filogenéticamente distantes —*Escherichia coli* (procariota) y *Schizosaccharomyces pombe* (eucariota)— permite distinguir mecanismos conservados del metabolismo central frente a respuestas específicas de cada linaje.

### Problema a resolver

Reproducir computacionalmente (y extender) el análisis de Vichi *et al.* mediante un **pipeline modular y validado** que, a partir de lecturas RNA-seq públicas, produzca:

- Matrices de conteos por gen.
- Resultados de expresión diferencial.
- Una comparación funcional entre organismos.

---

## Calendario de trabajo

> Fechas y plazos tentativos; pueden ajustarse durante el desarrollo del proyecto.

| Actividad                    | Fecha / periodo                       | Responsable                        | Entregable                              |
|------------------------------|---------------------------------------|------------------------------------|-----------------------------------------|
| Descripción de proyecto      | Septiembre (reuniones Ma y Vi 10–12) | Miryam Zamora, Alondra Márquez     | Documento `README.md`                   |
| Especificación de requisitos | Septiembre                            | Miryam y Alo                       | `README.md` con requisitos del pipeline |
| Análisis y diseño            | 20 de septiembre (aprox.)             | Miryam y Alo                       | Borrador de flujo / pseudocódigo       |
| Construcción                 | Octubre–noviembre                     | Miryam y Alo                       | Scripts en `src/`                       |
| Pruebas                      | Noviembre                             | Miryam y Alo                       | Informe breve de pruebas (markdown)     |
| Reporte de resultados        | Noviembre                             | Alo y Miryam                       | Documentos markdown + figuras           |
| Presentación del proyecto    | Diciembre                             | Alo y Miryam                       | Repositorio GitHub (release) + slides   |

---

## Preguntas de investigación

1. ¿Qué genes cambian significativamente entre respiración y fermentación en cada organismo (FDR < 0.05, \|log2FC\| ≥ 1)?
2. ¿Qué procesos/ontologías (p. ej., glucólisis, TCA, transporte de electrones, respuesta a estrés) están enriquecidos?
3. ¿Qué patrones son compartidos y cuáles son específicos entre *E. coli* y *S. pombe*?
4. ¿Cómo varían las métricas técnicas (tasa de mapeo, fracción asignada a genes, profundidad efectiva) entre condiciones y organismos?

---

## Hipótesis

- En **fermentación**, aumentarán genes relacionados con la glucólisis y disminuirán genes del ciclo de Krebs (TCA) y de la cadena de transporte de electrones (ETC); en **respiración**, se observará el patrón inverso.
- Existirá un **núcleo conservado** de respuestas metabólicas entre *E. coli* y *S. pombe*, pero con **diferencias regulatorias** propias de procariotas vs eucariotas (por ejemplo, organización operónica vs regulación por cromatina).

---

## Metodología

### Pasos generales

1. **Localización de la fuente de datos**  
   - Identificar BioProject/BioSample/GEO asociados al estudio (p. ej. PRJNA574477).
2. **Descarga de archivos de datos**  
   - Descarga de FASTQ desde SRA/ENA utilizando `sra-tools` o FTP/HTTPS.
3. **Inspección de datos**  
   - Verificar integridad de archivos.
   - Inspección básica (FastQC) de calidad de lecturas.
4. **Limpieza de datos**  
   - Recorte de adaptadores y baja calidad (Trim Galore!/cutadapt).
5. **Descripción de los datos**  
   - Estadísticos de mapeo, conteos por gen, distribución de conteos, PCA preliminar.
6. **Análisis de datos depurados**  
   - Generación de matrices de conteos.
   - Ejecución de DESeq2 (uno por organismo).
   - Análisis de enriquecimiento funcional (opcional).
7. **Obtención de resultados**  
   - Listas de DEGs, figuras (volcano, MA-plot, heatmaps, PCA).
   - Resumen comparativo *E. coli* vs *S. pombe*.

---

## Resultados esperados

Estructura esperada (conceptual) de salidas principales:

```text
data/metadata/samples.tsv

data/reference/*/{genome.fa, genes.gtf, index/}

results/align/<organism>/<sample>.sorted.{bam,bai}

results/counts/{Ecoli_counts.tsv, Spombe_counts.tsv}

results/DE/<organism>/{DEGs.tsv, normalized_counts.tsv}

results/plots/<organism>/{MA.png, volcano.png, PCA.png}

results/annotation/*.tsv

results/comparative/summary.tsv

logs/{download,validate,align,count,de}/*.log

docs/{presentacion_final.pdf, pipeline.svg}
```

--- 
## Análisis y diseño (pseudocódigo)

### Funciones principales

buscar_geo_por_proyecto(project_acc: str) -> dict
    # Busca el GEO asociado a un BioProject (p. ej., GSE)
    # Obtiene detalles para el primer GEO (o itera si hay >1)

obtener_sra_runinfo_por_bioproject(bioproject_acc: str) -> DataFrame
    # Usa esearch (db="sra") + efetch (rettype=runinfo)
    # Devuelve un DataFrame con columnas como:
    # Run, Experiment, Sample, BioSample, LibraryLayout,
    # scientific_name, ftp_path, etc.

seleccionar_fastq_urls(df_runs: DataFrame) -> list[str]
    # Extrae URLs FTP/HTTPS (o SRR IDs) para cada Run
    # Prioriza enlaces directos a *.fastq.gz o usa SRR para prefetch

buscar_ensamblajes_por_organismo(organism_term: str, filters: dict) -> DataFrame
    # Hace esearch en db="assembly" para el organismo dado
    # Aplica filtros (p. ej. refseq_only=True)

elegir_ensamblaje_referencia(df_assembly: DataFrame, preferencias: dict) -> row
    # Elige un ensamblado de referencia a partir de criterios:
    # - GCF preferido (GCF_000002945.2, GCF_000005845.2)
    # - Estado: "reference genome" o "representative genome"

listar_archivos_ftp(ftp_url: str) -> list[str]
    # Parsea la URL con urlparse
    # Convierte ftp://ftp.ncbi.nlm.nih.gov/... a https://...
    # Lista archivos disponibles en el directorio

seleccionar_y_descargar_archivos_genoma(row, outdir) -> local_paths
    # Usa FtpPath_RefSeq o FtpPath_GenBank del ensamblado seleccionado
    # Descarga genome.fa, genes.gff/gtf, etc., a outdir

descargar_fastq_con_sra_toolkit(run_ids: list[str], outdir: str) -> list[str]
    # Usa prefetch + fasterq-dump para obtener FASTQ
    # Comprime a fastq.gz y retorna rutas locales

alineamiento_por_replicado(fastq_pairs, reference_fasta, out_bam_dir) -> list[str]
    # Alinea cada par de FASTQ contra la referencia
    # Genera BAM ordenados + indexados en out_bam_dir

generar_matriz_de_conteos(bam_paths, annotation_gtf) -> "counts_matrix.tsv"
    # Cuenta lecturas por gen (featureCounts/HTSeq/coverageBed)
    # Combina resultados individuales en una matriz genes × muestras

analisis_DESeq2(counts_matrix, sample_metadata) -> "DE_results"
    # Construye objeto DESeq2, ajusta modelo, obtiene tabla de resultados
    # Aplica filtros por padj y |log2FC|

anotacion_y_comparacion_DEGs(DE_results_ecoli, DE_results_spombe) -> "reports"
    # Integra anotación funcional (GO/KEGG)
    # Compara patrones entre organismos y genera resúmenes/tablas/figuras
    
    # Flujo principal del pipeline
1. Obtener información GEO/SRA:
   - info_geo <- buscar_geo_por_proyecto("PRJNA574477")
   - df_runs  <- obtener_sra_runinfo_por_bioproject("PRJNA574477")
   - fastq_urls <- seleccionar_fastq_urls(df_runs)

2. Para cada organismo objetivo (S. pombe, E. coli):
   a) df_assembly <- buscar_ensamblajes_por_organismo(organism_term, filters={refseq_only: True})
   b) selected_assembly <- elegir_ensamblaje_referencia(df_assembly, preferencias={GCF_...})
   c) seleccionar_y_descargar_archivos_genoma(selected_assembly, BASE_DATA_DIR/{organism}/reference/)

3. Descargar FASTQ en BASE_DATA_DIR/fastq/ usando SRA toolkit o FTP/HTTPS.

4. Indexar referencias y alinear cada réplica:
   - alineamiento_por_replicado(fastq_pairs, reference_fasta, BASE_DATA_DIR/bam/)

5. Contar lecturas por gen:
   - generar_matriz_de_conteos(bam_paths, annotation_gtf) -> BASE_DATA_DIR/counts/

6. Ejecutar DESeq2 por organismo:
   - analisis_DESeq2(counts_matrix_ecoli, metadata_ecoli) -> BASE_DATA_DIR/results/Ecoli/
   - analisis_DESeq2(counts_matrix_spombe, metadata_spombe) -> BASE_DATA_DIR/results/Spombe/

7. Anotar y comparar resultados:
   - anotacion_y_comparacion_DEGs(DE_results_ecoli, DE_results_spombe) -> BASE_DATA_DIR/analysis/



