# Nombre del proyecto 
## High-throughput transcriptome sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth


**Fecha: 29/09/2025**
**Autores:** Miryam Zamora, Alondra Márquez

## Descripción del proyecto 

El presente proyecto tiene como propósito construir un pipeline reproducible de RNA-seq para ambos organismos y comparar firmas de expresión, así como y analizar, mediante herramientas bioinformáticas, los resultados reportados en el estudio:

>Vichi J., Salazar E., Jiménez Jacinto V., Olvera Rodríguez L., Grande R., Dantán-González E., Morett E., Hernández-Mendoza A. (2021). High-throughput transcriptome >sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth. PLOS ONE. https://doi.org/10.1371/journal.pone.0248513

En dicho trabajo, se llevó a cabo la secuenciación y análisis comparativo del transcriptoma de Escherichia coli (procariota) y Schizosaccharomyces pombe (eucariota unicelular) bajo condiciones de crecimiento respiratorias y fermentativas. A partir de estos datos, se identificaron genes diferencialmente expresados y se exploraron las similitudes y diferencias en los procesos metabólicos y regulatorios entre ambos organismos.

En este proyecto se emplearán directamente los datos proporcionados por los autores para llevar a cabo el pipeline de análisis bioinformático.

## Objetivo

Analizar los datos transcriptómicos de E. coli y S. pombe bajo condiciones de respiración y fermentación con el fin de identificar genes diferencialmente expresados, explorar su anotación funcional y realizar una comparación entre organismos que permita resaltar procesos conservados y divergentes.

    ## Objetivos específicos 

  1) Descargar los **12 FASTQ**   
  2) Descargar y versionar **genomas de referencia: Schizosaccharomyces pombe (fission yeast)     ASM294v3   GCF_000002945.2 reference genome; Escherichia coli str. K-12 substr. MG1655 (E. coli)     ASM584v2   GCF_000005845.2 reference genome** 
  (FASTA) y **anotaciones** (GTF/GFF).  (Que se descargan automáticamente a ../../data/)
  3) Alinear por réplica y generar **BAM ordenados + BAI**.  
  4) Obtener **matriz de conteos por gen**.  
  5) Realizar **DESeq2** por organismo.  
  6) **Anotar DEGs** y **comparar** patrones compartidos/divergentes.

## Planteamiento del problema

Contexto. La respiración y la fermentación inducen reprogramaciones transcriptómicas profundas. Comparar estas respuestas en dos organismos filogenéticamente distantes—Escherichia coli (procariota) y Schizosaccharomyces pombe (eucariota)—permite distinguir mecanismos conservados del metabolismo central frente a respuestas específicas de linaje.

Problema a resolver. Reproducir computacionalmente (y extender) el análisis de Vichi et al. con un pipeline modular y validado que, a partir de lecturas RNA-seq públicas, produzca matrices de conteo, resultados de expresión diferencial y una comparación funcional entre organismos.

## Calendario de trabajo

[Definir de manera general la actividades que se requerirán para el proyecto. Por ejemplo:]

| Actividad | Fecha   | Responsable  | Entregable |
|----------|----------|----------|----------|
| Descripción de proyecto    | Reunirnos Martes y viernes de 10 a 12  | Miryam Zamora, Alondra Márquez  | Documento markdown-README.md |
| Especificación de requisitos    | septiembre   | Alo y Miryam   | Documento markdown-README.md   |
| Análisis y diseño   | 20 septiembre  |
| Construcción   | octubre, noviembre  |  Miryam y Alo    | Scripts |
| Pruebas   | noviembre  |  Miryam y Alo   | Documento markdown |
| Reporte de resultados  | noviembre  |  Alo y Miryam   | Documentos markdown |
| Presentación del proyecto   | diciembre  |  Alo y Miryam   | repositorio GitHub (release)|


**Preguntas de investigación**

¿Qué genes cambian significativamente entre respiración y fermentación en cada organismo (FDR < 0.05, |log2FC| ≥ 1)?

¿Qué procesos/ontologías (p. ej., glucólisis, TCA, transporte de electrones, estrés) están enriquecidos?

¿Qué patrones son compartidos y cuáles son específicos entre E. coli y S. pombe?

¿Cómo varían métricas técnicas (tasa de mapeo, fracción asignada a genes) entre condiciones y organismos?

Hipótesis.

En fermentación aumentarán genes glucolíticos y disminuirán genes del TCA/ETC; en respiración ocurrirá el patrón inverso.

Existirá un núcleo conservado de respuestas metabólicas, con diferencias regulatorias propias de procariotas vs eucariotas.


## Metodología ## 

Paso 1: Localización de fuente de datos  
Paso 2: Descarga de archivos de datos  
Paso 3: Inspección de datos  
Paso 4: Limpieza de datos
Paso 5: Descripción de los datos  
Paso 5: Análsis de datos depurados  
Paso 6: Obtención de resultados  



## Resultados esperados

data/metadata/samples.tsv

data/reference/*/{genome.fa, genes.gtf, index/}

results/align/*/*.sorted.{bam,bai}

results/counts/{Ecoli_counts.tsv, Spombe_counts.tsv}

results/DE/*/{DEGs.tsv, normalized_counts.tsv}

results/plots/*/{MA.png, volcano.png, PCA.png}

results/annotation/*.tsv

results/comparative/summary.tsv

logs/{download,validate,align,count,de}/*.log

docs/{presentacion_final.pdf, pipeline.svg}



## Análisis y Diseño

Funciones (pseudocódigo)
1. buscar_geo_por_proyecto(project_acc: str) -> dict
    #Busca GEO asociado al BioProject (p. ej. GSE) 
    #Obtener detalles para el primer GEO (o iterar si hay >1)

2. obtener_sra_runinfo_por_bioproject(bioproject_acc: str) -> DataFrame 
    #Usa esearch (sra) + efetch rettype=runinfo para obtener runinfo CSV
    # Campos útiles: Run, Experiment, Sample, BioSample, LibraryLayout, scientific_name, ftp_path

3. seleccionar_fastq_urls(df_runs: DataFrame) -> list[str]
    #Extrae FTP https URLs (o SRA identifiers) para cada Run
    #Priorizar: ftp://.../fastq.gz (si existe) o usar run accession para prefetch

4. buscar_ensamblajes_por_organismo(organism_term: str, filters: optional) -> DataFrame
    #Hacer esearch por organimso en db="assembly"

5. elegir_ensamblaje_referencia(df_assembly: DataFrame, preferencias: dict) -> row

6. listar_archivos_ftp(ftp_url: str) -> list[str]
    # Parsear url con urlparse
    # Convertir ftp://ftp.ncbi.nlm.nih.gov/... a https

7. seleccionar_y_descargar_archivos_genoma(row, outdir) -> local_paths
# row debe contener FtpPath_RefSeq o FtpPath_GenBank

8. descargar_fastq_con_sra_toolkit(run_ids: list[str], outdir) -> local_fastq_paths

9. alineamiento_por_replicado(fastq_pairs, reference_fasta, out_bam_dir) -> bam_paths

10. generar_matriz_de_conteos(bam_paths, annotation_gtf) -> counts_matrix.tsv

11. análisis_DESeq2(counts_matrix, sample_metadata) -> DE_results

12. anotación_y_comparación_DEGs(DE_results_ecoli, DE_results_spombe) -> reports


- **Flujo principal**:

	1. Obtener GEO/SRA info:
   - info_geo <- buscar_geo_por_proyecto("PRJNA574477" o ProjectAcc)
   - df_runs <- obtener_sra_runinfo_por_bioproject("PRJNA574477")

2. Obtener lista de FASTQ (por run) -> seleccionar_fastq_urls(df_runs)

3. Para cada organismo objetivo (S. pombe, E. coli):
   a. df_assembly <- buscar_ensamblajes_por_organismo(organism_term, filters.refseq_only=True)
   b. selected_assembly <- elegir_ensamblaje_referencia(df_assembly, preferencia por GCF_000002945.2 / GCF_000005845.2 si se especifica)
   c. descargar archivos genoma y anotación -> guardar en BASE_DATA_DIR/{organism}/reference/

4. Descargar FASTQ en BASE_DATA_DIR/fastq/ (usar SRA toolkit o FTP HTTPS)
5. Indexar referencias y alinear cada réplica -> generar BAM + BAI en BASE_DATA_DIR/bam/
6. Contar lecturas por gen -> BASE_DATA_DIR/counts/
7. Ejecutar DESeq2 por organismo -> resultados en BASE_DATA_DIR/results/
8. Anotar y comparar resultados -> informes y visualizaciones en BASE_DATA_DIR/analysis/
9. Guardar metadata, versiones (AssemblyAccession, Fecha, Hash archivos) en BASE_DATA_DIR/metadata/

## DIAGRAMA - DESCRIPCIÓN DEL PIPELINE

[## Flujo general del pipeline

**Entrada:** Datos públicos RNA-seq (SRA/ENA)  

[ Datos públicos RNA-seq (SRA/ENA) ]
↓
[ Diseño de búsqueda ]
(egquery / espell / einfo)
↓
[ IDs SRA ] ← esearch (bioproject / sra)
↓
[ RunInfo.csv + samples.tsv ] ← efetch (runinfo)
↓
[ Descarga FASTQ ]

+----------------- Referencias -----------------+
| esearch / esummary (assembly) |
| elink (assembly → nuccore) |
| efetch (fasta / gff) |
+-----------------------------------------------+
↓
[ genome.fa + genes.gtf ]
        ↓
[ Alineamiento (HISAT2 / STAR) ]
↓
[ Conteos (featureCounts / HTSeq) ]
↓
[ DESeq2 por organismo ]
↓
[ Anotación de DEGs ]
↓
[ Comparación: E. coli vs S. pombe ]
