# High-throughput transcriptome sequencing and comparative analysis of *Escherichia coli* and *Schizosaccharomyces pombe* in respiratory and fermentative growth

**Fecha:** 06/10/2025  
**Autores:** Miryam Zamora, Alondra Márquez

----

## Descripción

Reproducimos y extendemos el análisis de Vichi _et al._ (2021) para comparar las respuestas transcriptómicas de _E. coli_(procariota) y _S. pombe_ (eucariota) en respiración vs fermentación.

Se construyó un **pipeline reproducible de RNA-seq**, implementado en Python y ejecutable en clúster (SGE), que:

- Descarga datos públicos (SRA/ENA) y genomas de referencia.
- Organiza las lecturas por organismo y muestra.
- Limpia y recorta lecturas (Trim Galore!).
- Alinea lecturas con **BWA**
- Cuantifica expresión génica por cobertura (bedtools `coverage` / `coverageBed`).
- Construye tablas de cobertura genes × muestras.
- Ejecuta análisis de expresión diferencial (DESeq2/pydeseq2) por organismo.
- Clasifica genes en **UP/DOWN** regulados.
- Genera gráficas (volcano, barplots, heatmaps, PCA, etc.).
- Organiza resultados en una estructura de carpetas clara y trazable.

---
## Objetivos

### Objetivo general

Identificar genes diferencialmente expresados y patrones conservados/divergentes entre _E. coli_ y _S. pombe_ bajo respiración y fermentación, usando un flujo reproducible de RNA-seq.

### Objetivos específicos

- Descargar las 12 bibliotecas FASTQ (2 organismos × 2 condiciones × 3 réplicas).
    
- Descargar y versionar genoma (+ anotación) para:
    
    - _E. coli_ — `GCF_000005845.2`
        
    - _S. pombe_ — `GCF_000002945.2`
        
- Alineación por réplica usando BWA-MEM.
    
- Generar tablas de cobertura por gen (genes × muestras).
    
- Realizar análisis DESeq2 por organismo.
    
- Clasificar genes diferencialmente expresados (UP/DOWN).
    
- Generar figuras de apoyo (QC, volcano, heatmaps, PCA).
    
- Comparar resultados entre organismos.

----

## Descripción del pipeline

### Diagrama ASCII (vista general)

```
[Descarga de metadatos y genomas]
       |
       v
[Descarga de FASTQ crudos desde SRA/ENA]
       |
       v
[Limpieza y recorte de FASTQ (Trim Galore!)]
       |
       v
[Alineamiento BWA-MEM → BAM ordenados + indexados]
       |
       v
[Cuantificación por gen (bedtools coverage) → *.count.txt]
       |
       v
[Construcción tabla cobertura genes × muestras (Ecoli/Spombe)]
       |
       v
[DESeq2 / pydeseq2 → resultados crudos de DE]
       |
       v
[Clasificación de DEGs (UP/DOWN) + listas genes_UP/genes_DOWN]
       |
       v
[Gráficas: volcano, barplots, heatmaps, PCA, boxplots, densidades]

```

---

## Estructura de carpetas (vista real del proyecto)

A partir de la estructura actual, las carpetas clave son:

```
.
├── data
│   ├── BWA
│   │   ├── e_coli_bam/          # BAM + índices + *.count.txt de E. coli
│   │   └── s_pombe_bam/         # BAM + índices + *.count.txt de S. pombe
│   ├── Escherichia_coli_GCF_000005845.2/      # Genoma + GFF + CDS + rRNA
│   ├── Schizosaccharomyces_pombe_GCF_000002945.2/ # Genoma + GFF + CDS + rRNA
│   ├── Escherichia_coli_str._K-12_substr._MG1655/ # FASTQ crudos E. coli
│   ├── Schizosaccharomyces_pombe/                  # FASTQ crudos S. pombe
│   ├── metadata/
│   │   └── samples.csv           # metadatos generales
│   ├── E.coli_samples.csv        # diseño específico E. coli
│   ├── S.pombe_samples.csv       # diseño específico S. pombe
│   ├── metadata_ecoli_defseq.csv # diseño para DESeq2 (E. coli)
│   └── metadata_spombe_defseq.csv# diseño para DESeq2 (S. pombe)
├── results
│   ├── Cleaned_fastq/            # FASTQ recortados (Trim Galore!)
│   │   ├── Escherichia_coli_str._K-12_substr._MG1655/
│   │   └── Schizosaccharomyces_pombe/
│   ├── analisis_preliminar/      # boxplots, densidades y PCA por organismo
│   ├── DE/
│   │   ├── DESeq2_Ecoli/         # resultados DESeq2 E. coli
│   │   └── DESeq2_Spombe/        # resultados DESeq2 S. pombe
│   ├── figuras/
│   │   ├── E_coli/               # volcano, barplots, heatmaps E. coli
│   │   └── S_pombe/              # volcano, barplots, heatmaps S. pombe
│   ├── tabla_cobertura_Ecoli.csv # matriz genes × muestras (E. coli)
│   └── tabla_cobertura_Spombe.csv# matriz genes × muestras (S. pombe)
├── logs/                         # archivos de log del pipeline
├── src/
│   ├── main_download.py          # descarga referencias + FASTQ
│   ├── clean_fastq.py            # recorte de FASTQ (Trim Galore!)
│   ├── main_aligments_counts.py  # BWA + cobertura + tablas de cobertura
│   ├── analisis_preliminar.py    # boxplot, densidad, PCA
│   ├── main_deseq_fit.py         # ajuste de modelo DESeq2/pydeseq2
│   ├── main_deseq_clasificar.py  # clasificación de DEGs UP/DOWN
│   ├── main_plots.py             # volcano, barplots, heatmaps
│   ├── stats.py                  # utilidades estadísticas
│   └── modules/
│       ├── bioproject.py         # manejo de BioProjects/Entrez
│       ├── sra_download.py       # descarga de SRR/FASTQ
│       ├── genome_reference.py   # descarga + manejo de genomas
│       ├── cobertura.py          # BWA + bedtools coverage
│       ├── expresion_diff.py     # funciones para DESeq2/pydeseq2
│       └── plots_DE.py           # funciones para gráficas de DE
└── docs/                         # documentos extra (presentación, etc.)

```

### ¿Qué hay en cada parte?

- `data/Escherichia_coli_GCF_000005845.2/` y  
    `data/Schizosaccharomyces_pombe_GCF_000002945.2/`  
    → Genomas de referencia, anotaciones GFF y archivos derivados (índices para BWA).
    
- `data/Escherichia_coli_str._K-12_substr._MG1655/` y  
    `data/Schizosaccharomyces_pombe/`  
    → FASTQ **crudos** descargados desde SRA/ENA, organizados por **GSM**/SRR.
    
- `results/Cleaned_fastq/`  
    → FASTQ **recortados y filtrados** por Trim Galore! (`*_val_1.fq.gz`, `*_val_2.fq.gz`)
    
    - reportes de trimming (`*_fastq.gz_trimming_report.txt`).
        
- `data/BWA/e_coli_bam/` y `data/BWA/s_pombe_bam/`  
    → Alineamientos BAM (`*.fq.sorted.bam`), índices `.bai` y archivos `*.count.txt` generados por `bedtools coverage`.
    
- `results/tabla_cobertura_Ecoli.csv` y `results/tabla_cobertura_Spombe.csv`  
    → Tablas de cobertura **genes × muestras** que alimentan el análisis DE.
    
- `results/DE/DESeq2_Ecoli` y `results/DE/DESeq2_Spombe`  
    → Resultados de expresión diferencial:
    
    - `deseq2_resultados_crudos.csv`
        
    - `counts_filtrados_cpm.csv`
        
    - Carpetas `*_clasificado/` con:
        
        - `deseq2_resultados_clasificados.csv`
            
        - `genes_UP.txt`
            
        - `genes_DOWN.txt`
            
- `results/figuras/`  
    → Volcano plots, barplots y heatmaps por organismo.
    
- `results/analisis_preliminar/`  
    → QC de conteos: boxplots, densidades, PCA y resúmenes por organismo.
    
- `logs/`  
    → Trazas de ejecución (`pipeline_master.log`, `*_nohup.out`, logs de DESeq, cobertura, trimming…).

---

## Flujo general: Entradas → Scripts → Salidas (con rutas reales)

### 1. Descarga de referencias y FASTQ (`src/main_download.py`)

- **Entrada:** IDs de BioProject/SRR definidos en el script.
    
- **Salida principal:**
    
    - Genomas y anotaciones en:
        
        - `data/Escherichia_coli_GCF_000005845.2/`
        - `data/Schizosaccharomyces_pombe_GCF_000002945.2/`
            
    - FASTQ crudos organizados en:
        
        - `data/Escherichia_coli_str._K-12_substr._MG1655/`
        - `data/Schizosaccharomyces_pombe/`
            
    - Metadatos:
        
        - `data/metadata/samples.csv`
        - `data/E.coli_samples.csv`
        - `data/S.pombe_samples.csv`

### 2. Limpieza de FASTQ (`src/clean_fastq.py`)

- **Entrada:**
    
    - FASTQ crudos (`data/Escherichia_coli_str._K-12_substr._MG1655/*_1.fastq.gz`, etc.)
- **Herramientas:** Trim Galore! (basado en `cutadapt` / `fastqc`).
- **Salida:**
    - FASTQ recortados en:
        - `results/Cleaned_fastq/Escherichia_coli_str._K-12_substr._MG1655/*_val_*.fq.gz`
        - `results/Cleaned_fastq/Schizosaccharomyces_pombe/*_val_*.fq.gz`
    - Reportes de trimming (`*_trimming_report.txt`).

### 3. Alineamiento + cobertura (`src/main_aligments_counts.py` + `modules/cobertura.py`)

- **Entrada:**
    
    - FASTQ limpios (`results/Cleaned_fastq/...`)
    - Genomas e índices BWA en:
        - `data/Escherichia_coli_GCF_000005845.2/`
        - `data/Schizosaccharomyces_pombe_GCF_000002945.2/`
            
- **Herramientas:**
    
    - `bwa mem`
    - `samtools sort` / `samtools index`
    - `bedtools coverage` (`coverageBed`)
        
- **Salida:**
    
    - E. coli:
        
        - `data/BWA/e_coli_bam/*.fq.sorted.bam`
        - `data/BWA/e_coli_bam/*.fq.sorted.bam.bai`
        - `data/BWA/e_coli_bam/*.count.txt`
            
    - S. pombe:
        
        - `data/BWA/s_pombe_bam/*.fq.sorted.bam`
        - `data/BWA/s_pombe_bam/*.fq.sorted.bam.bai`
        - `data/BWA/s_pombe_bam/*.count.txt`
            
    - Tablas combinadas genes × muestras:
        
        - `results/tabla_cobertura_Ecoli.csv`
        - `results/tabla_cobertura_Spombe.csv`

### 4. Análisis preliminar de conteos (`src/analisis_preliminar.py`)

- **Entrada:**
    
    - `results/tabla_cobertura_Ecoli.csv`
    - `results/tabla_cobertura_Spombe.csv`
        
- **Salida:**
    
    - `results/analisis_preliminar/Ecoli/{boxplot_Ecoli.png, densidad_Ecoli.png, pca_Ecoli.png, resumen_Ecoli.csv}`
    - `results/analisis_preliminar/Spombe/{boxplot_Spombe.png, densidad_Spombe.png, pca_Spombe.png, resumen_Spombe.csv}`

### 5. Ajuste del modelo de expresión diferencial (`src/main_deseq_fit.py`)

- **Entrada:**
    
    - Tablas de cobertura:
        - `results/tabla_cobertura_Ecoli.csv`
        - `results/tabla_cobertura_Spombe.csv`
            
    - Diseños experimentales:
        - `data/metadata_ecoli_defseq.csv`
        - `data/metadata_spombe_defseq.csv`
            
- **Herramientas:**
    - `DESeq2` (vía R) o `pydeseq2` (Python), según implementación en `modules/expresion_diff.py`.
        
- **Salida:**
    
    - `results/DE/DESeq2_Ecoli/deseq2_resultados_crudos.csv`
    - `results/DE/DESeq2_Ecoli/counts_filtrados_cpm.csv`
    - `results/DE/DESeq2_Spombe/deseq2_resultados_crudos.csv`
    - `results/DE/DESeq2_Spombe/counts_filtrados_cpm.csv`

### 6. Clasificación de DEGs (`src/main_deseq_clasificar.py`)

- **Entrada:**
    
    - `results/DE/DESeq2_*/deseq2_resultados_crudos.csv`
- **Salida:**
    
    - E. coli:
        
        - `results/DE/DESeq2_Ecoli/DESeq2_Ecoli_clasificado/`
            - `deseq2_resultados_clasificados.csv`
            - `genes_UP.txt`
            - `genes_DOWN.txt`
                
    - S. pombe:
        
        - `results/DE/DESeq2_Spombe/DESeq2_Spombe_clasificado/`
            - `deseq2_resultados_clasificados.csv`
            - `genes_UP.txt`
            - `genes_DOWN.txt`

### 7. Gráficas de DEGs (`src/main_plots.py` + `modules/plots_DE.py`)

- **Entrada:**
    
    - Resultados clasificados y tablas de conteos filtrados.
- **Salida:**
    - `results/figuras/E_coli/{volcano_E._coli.png, barplot_DE_E._coli.png, heatmap_top50_E._coli.png, heatmap_top100_E._coli.png}`
    - `results/figuras/S_pombe/{volcano_S.pombe.png, barplot_DE_S.pombe.png, heatmap_top50_S.pombe.png}`

---

## Metadatos de muestras

Los metadatos generales se encuentran en:  
`data/metadata/samples.csv`

Además, se utilizan diseños específicos para DESeq2:

- `data/metadata_ecoli_defseq.csv`
- `data/metadata_spombe_defseq.csv`

| sample_id   | organism | condition | replicate | read_type | GSM        | SRR         | genome_ref       | control |
|-------------|----------|-----------|-----------|-----------|------------|-------------|------------------|---------|
| GSM4099077  | Ecoli    | Resp      | 1         | PE        | GSM4099077 | SRR10192862 | GCF_000005845.2  | yes     |
| GSM4099078  | Ecoli    | Resp      | 2         | PE        | GSM4099078 | SRR10192863 | GCF_000005845.2  | yes     |
| GSM4099079  | Ecoli    | Resp      | 3         | PE        | GSM4099079 | SRR10192864 | GCF_000005845.2  | yes     |
| ...         | ...      | ...       | ...       | ...       | ...        | ...         | ...              | ...     |


---

## Requisitos de software

| Categoría              | Herramientas / Versiones sugeridas                                              |
| ---------------------- | ------------------------------------------------------------------------------- |
| **Python (3.10+)**     | `biopython`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pydeseq2` (si se usa) |
| **Descarga SRA**       | `SRA Toolkit` (`prefetch`, `fasterq-dump`)                                      |
| **Trimming/QC**        | `Trim Galore!`, `cutadapt`, `FastQC`, `MultiQC`                                 |
| **Alineamiento**       | `bwa`, `samtools`                                                               |
| **Cobertura**          | `bedtools` (`coverage`)                                                         |
| **R (4.2+)**           | `DESeq2`, `ggplot2`, `pheatmap`, `tidyverse`                                    |
| **Cluster (opcional)** | SGE (`qsub`) u otro scheduler                                                   |

---

## Cómo reproducir el análisis (versión scripts)

Asumiendo que estás en la raíz del repositorio:

### 0. Clonar el repositorio

```bash 
git clone https://github.com/USER/REPO.git
cd REPO
```

### 1. Crear entorno
```bash 
mamba env create -f environment.yml    # si existe mamba activate ecoli_spombe_rnaseq`
```

### 2. Descarga de genomas y FASTQ

```bash
python src/main_download.py
```

- Descarga genomas a `data/Escherichia_coli_GCF_000005845.2/` y `data/Schizosaccharomyces_pombe_GCF_000002945.2/`.
- Descarga FASTQ crudos a las carpetas de organismo correspondientes.
- Genera/actualiza metadatos (`data/metadata/samples.csv`, etc.).

### 3. Limpieza de FASTQ

```bash
python src/clean_fastq.py
```


- Aplica Trim Galore! sobre los FASTQ crudos.
- Escribe FASTQ limpios en `results/Cleaned_fastq/...`.

### 4. Alineamiento y conteos por gen

```
python src/main_aligments_counts.py
```


- Alinea con BWA-MEM, ordena e indexa BAM.
- Ejecuta `bedtools coverage` para obtener `*.count.txt`.
- Combina conteos en:
    - `results/tabla_cobertura_Ecoli.csv`
    - `results/tabla_cobertura_Spombe.csv`


### 5. Análisis preliminar de conteos (opcional pero recomendado)

```
python src/analisis_preliminar.py
```

- Genera boxplots, densidades, PCA y resúmenes en `results/analisis_preliminar/`.


### 6. Ajuste del modelo DESeq2/pydeseq2

```

`python src/main_deseq_fit.py \   --counts-ecoli results/tabla_cobertura_Ecoli.csv \   --counts-spombe results/tabla_cobertura_Spombe.csv \   --metadata-ecoli data/metadata_ecoli_defseq.csv \   --metadata-spombe data/metadata_spombe_defseq.csv`
```

- Guarda resultados crudos y conteos filtrados en `results/DE/DESeq2_*`.

### 7. Clasificación de DEGs

```
python src/main_deseq_clasificar.py
```


- Genera `genes_UP.txt`, `genes_DOWN.txt` y `deseq2_resultados_clasificados.csv` para cada organismo en `results/DE/DESeq2_*/..._clasificado/`.

### 8. Gráficas finales de DE

`python src/main_plots.py`

- Produce volcano plots, barplots y heatmaps en `results/figuras/E_coli/` y `results/figuras/S_pombe/`.

---
## Buenas prácticas y reproducibilidad

- **Logs:** revisar `logs/*.log` para diagnosticar problemas (trimming, cobertura, DESeq2, etc.).
- **Versionado de datos:** anotar accesiones SRA y versiones de genomas en un archivo `docs/VERSIONS.txt`. 
- **Seeds:** fijar semillas en DESeq2/pydeseq2 para mejorar reproducibilidad de ciertos pasos (ej. PCA).

---
## Resultados esperados (resumen)

- Tablas de DEGs por organismo (con log2FoldChange, padj, etc.).
- Listas de genes UP/DOWN:
    - `results/DE/DESeq2_*/.../genes_UP.txt`
        
    - `results/DE/DESeq2_*/.../genes_DOWN.txt`
        
- Gráficas:
    - QC de conteos (boxplots, densidad, PCA).
        
    - Volcano plots, barplots y heatmaps de genes más significativos.
        
- Tablas de cobertura por gen:
    - `results/tabla_cobertura_Ecoli.csv`
        
    - `results/tabla_cobertura_Spombe.csv`

---
## Licencia


---
## Contacto

- **Miryam Zamora** — correo institucional
    
- **Alondra Márquez** — alondram@lcg.unam.mx

Para dudas o sugerencias relacionadas con el código o el análisis, abre un _issue_ en el repositorio o contacta por correo.
---

Cita
Vichi J., Salazar E., Jiménez Jacinto V., Olvera Rodríguez L., Grande R., Dantán-González E., Morett E., Hernández-Mendoza A. (2021).
High-throughput transcriptome sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth.
PLOS ONE. doi:10.1371/journal.pone.0248513
