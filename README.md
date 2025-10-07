High-throughput transcriptome sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth

Fecha: 06/10/2025
Autores: Miryam Zamora, Alondra Márquez

Descripción

Reproducimos y extendemos el análisis de Vichi et al. (2021) para comparar las respuestas transcriptómicas de E. coli (procariota) y S. pombe (eucariota) en respiración vs fermentación. Construimos un pipeline reproducible que descarga datos públicos, alinea, cuantifica, realiza DESeq2 y compara DEGs entre organismos.

Objetivo general

Identificar genes diferencialmente expresados y patrones conservados/divergentes entre E. coli y S. pombe bajo respiración y fermentación, usando un flujo reproducible de RNA-seq.

Objetivos específicos

Descargar las 12 bibliotecas FASTQ (2 org × 2 condiciones × 3 réplicas).

Descargar y versionar genoma (+ anotación) para E. coli (GCF_000005845.2) y S. pombe (GCF_000002945.2).

Alinear cada réplica por condición y organismo.

Generar una matriz de conteos por gen.

Realizar DESeq2 por organismo.

Anotar los DEGs.

Comparar resultados: similitudes/diferencias entre organismos.

Descripción del pipeline

Diagrama (ASCII)

[SRA/ENA FASTQ (12)]
        |
        v
[Descarga FASTQ]
        |
        v                      +---------------- Referencias ---------------+
[Alineamiento (HISAT2/STAR)]   | assemblies -> (nuccore) -> FASTA + GFF     |
        |                      | genome.fa + genes.gff                      |
        |                      +--------------------------------------------+
        v
[Conteo por gen (featureCounts/HTSeq)]
        |
        v
[DESeq2 por organismo]
        |
        v
[Anotación de DEGs]
        |
        v
[Comparación E. coli vs S. pombe]


Entradas → Herramientas → Salidas

Referencias (genoma + anotación)

Herramientas: Biopython/Entrez (esearch, esummary), descarga FTP/HTTPS.

Salida:

data/E_coli_MG1655/{genome.fa, genes.gff}

data/S_pombe_ASM294v3/{genome.fa, genes.gff}

VERSIONS.txt con accesiones y checksums.

Lecturas (12 FASTQ)

Herramientas: SRA Toolkit (prefetch, fasterq-dump) o descarga desde ENA.

Salida: data/fastq/raw/<sample>_R{1,2}.fastq.gz.

Alineamiento

Herramientas: HISAT2 (o STAR) + Samtools.

Salida: results/align/<org>/<sample>.sorted.bam + .bai.

Conteos

Herramientas: featureCounts (Subread) o HTSeq.

Salida: results/counts/{Ecoli_counts.tsv, Spombe_counts.tsv}.

DE

Herramientas: R/DESeq2.

Salida: results/DE/<org>/{DEGs.tsv, normalized_counts.tsv} + gráficos (MA.png, volcano.png, PCA.png).

Anotación

Herramientas: Entrez/archivos de referencia.

Salida: results/annotation/{Ecoli_annotation.tsv, Spombe_annotation.tsv}.

Comparación

Herramientas: scripts propios (Python/R).

Salida: results/comparative/summary.tsv (intersecciones, pathways comunes, etc.).

Requisitos (software)

Python 3.10+: biopython, pandas, numpy

SRA Toolkit: prefetch, fasterq-dump

Alineador: HISAT2 (o STAR)

Samtools

Subread (featureCounts) o HTSeq

R 4.2+: DESeq2, ggplot2

(Opcional) gffread para convertir GFF→GTF, FastQC/MultiQC para QC

Estructura de carpetas
.
├                    
├── src/                           # Donde ejecutaremos los main
│   ├──modules                     # scripts con funciones
├── data/                           
│   ├── fastq/
│   │   └── raw/                   # *.fastq.gz
│   ├── E_coli_MG1655/             # genome.fa, genes.gff
│   └── S_pombe_ASM294v3/          # genome.fa, genes.gff
├── results/
│   ├── align/
│   ├── counts/
│   ├── DE/                         #Salidas de DESeq2
│   ├── annotation/
│   └── comparative/
├── docs/
│   ├── pipeline_ascii.md
│   └── presentacion_final.pdf
├── test/                          #códigos de prueba
└── README.md

Metadatos de muestras

Guarda esto como data/metadata/samples.tsv y sustituye las SRR reales. Si tu librería es single-end, cambia read_type a SR.

sample_id	organism	condition	replicate	read_type	accession	bioproject	genome_ref	strand	control	tissue
Ecoli_R1	Ecoli	Resp	1	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000005845.2	0	yes	cultivo; células completas
Ecoli_R2	Ecoli	Resp	2	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000005845.2	0	yes	cultivo; células completas
Ecoli_R3	Ecoli	Resp	3	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000005845.2	0	yes	cultivo; células completas
Ecoli_F1	Ecoli	Ferm	1	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000005845.2	0	no	cultivo; células completas
Ecoli_F2	Ecoli	Ferm	2	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000005845.2	0	no	cultivo; células completas
Ecoli_F3	Ecoli	Ferm	3	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000005845.2	0	no	cultivo; células completas
Spombe_R1	Spombe	Resp	1	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000002945.2	0	yes	cultivo; células completas
Spombe_R2	Spombe	Resp	2	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000002945.2	0	yes	cultivo; células completas
Spombe_R3	Spombe	Resp	3	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000002945.2	0	yes	cultivo; células completas
Spombe_F1	Spombe	Ferm	1	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000002945.2	0	no	cultivo; células completas
Spombe_F2	Spombe	Ferm	2	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000002945.2	0	no	cultivo; células completas
Spombe_F3	Spombe	Ferm	3	PE	SRRxxxxx	PRJNAxxxxxx	GCF_000002945.2	0	no	cultivo; células completas

Cómo ejecutar (resumen)

Partimos de las descargas del genoma de referencias a data/E_coli_MG1655 y data/S_pombe_ASM294v3.

1) Descarga FASTQ (ejemplos SRA Toolkit)

# Por SRR (prefetch + fasterq-dump → gzip)
prefetch SRRxxxxxx
fasterq-dump SRRxxxxxx -O data/fastq/raw -e 8
gzip data/fastq/raw/*.fastq


2) Indexar referencias (HISAT2)

hisat2-build data/E_coli_MG1655/genome.fa data/E_coli_MG1655/index/genome
hisat2-build data/S_pombe_ASM294v3/genome.fa data/S_pombe_ASM294v3/index/genome


3) Alinear (ejemplo PE)

hisat2 -x data/E_coli_MG1655/index/genome -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  | samtools sort -o results/align/Ecoli/sample.sorted.bam
samtools index results/align/Ecoli/sample.sorted.bam


4) Conteos (featureCounts)

featureCounts -T 8 -a data/E_coli_MG1655/genes.gff -o results/counts/Ecoli_counts.tsv \
  -t gene -g ID results/align/Ecoli/*.sorted.bam

featureCounts -T 8 -a data/S_pombe_ASM294v3/genes.gff -o results/counts/Spombe_counts.tsv \
  -t exon -g gene_id results/align/Spombe/*.sorted.bam


5) DESeq2 (R)

# scripts/deseq2.R (ejecutar con Rscript)
# Lee counts + diseño, corre DESeq2 por organismo, escribe DEGs y gráficas.


6) Anotación y comparación

Python/R: fusionar DEGs con anotación (GFF) y producir results/comparative/summary.tsv.

Validación (como pide la rúbrica)

Verificación de accesiones y versiones: VERSIONS.txt en cada referencia.

Confirmar layout (PE/SR) y strandness (RSeQC/featureCounts summary).

Controles: Resp = control; Ferm = tratamiento.

Consistencia de réplicas: 3 por grupo; matriz de conteos sin NA.

Licencia

Este repositorio se distribuye bajo la licencia MIT (ver LICENSE). Los datos originales pertenecen a sus respectivos autores y repositorios.

Cita

Vichi J., Salazar E., Jiménez Jacinto V., Olvera Rodríguez L., Grande R., Dantán-González E., Morett E., Hernández-Mendoza A. (2021). High-throughput transcriptome sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth. PLOS ONE. doi:10.1371/journal.pone.0248513