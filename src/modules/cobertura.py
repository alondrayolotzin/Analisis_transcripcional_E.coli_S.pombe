"""
cobertura.py
Autor: Miryam Zamora
Fecha: 2025-11-23

Funciones para:
  - Alinear lecturas FASTQ pareadas con BWA MEM y generar BAM ordenados e indexados.
  - Cuantificar lecturas por gen/región anotada usando bedtools coverage (coverageBed).
  - Leer múltiples archivos de conteo (*.count.txt), normalizarlos y combinarlos
    en una tabla de cobertura genes x muestras lista para usar en análisis posteriores
    (por ejemplo, DESeq2).

Este módulo es llamado por main.cobertura.py y forma parte del pipeline de RNA-seq.

=================
"""
import pandas as pd
import os
import re
import subprocess
from pathlib import Path
import gzip
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

def alinear_con_bwa(bwa_index, reads_dir, out_dir, threads=8, sort_mem="2G", tmpdir="/dev/shm"):
    """
    
    Alinea lecturas pareadas (R1/R2) con BWA MEM y produce archivos BAM
    ordenados e indexados con samtools.

    Flujo general
    -------------
    1. Se asegura de que el genoma de referencia esté en formato FASTA sin comprimir.
       - Si el archivo termina en .gz, se descomprime temporalmente.
    2. Verifica que el índice de BWA exista para ese FASTA; si no, lo genera.
    3. Busca archivos FASTQ comprimidos (*.fastq.gz, *.fq.gz) en el directorio indicado.
    4. Agrupa los archivos en pares R1/R2 utilizando patrones en el nombre.
    5. Para cada par:
       - Ejecuta `bwa mem` para alinear.
       - Pasa la salida a `samtools view` para convertir a BAM.
       - Ordena el BAM con `samtools sort`.
       - Indexa el BAM final con `samtools index`.

    Parámetros
    ----------
    bwa_index : str o pathlib.Path
        Ruta al genoma de referencia. Puede ser un FASTA (.fa/.fna)
        o un FASTA comprimido (.fa.gz / .fna.gz).
    reads_dir : str o pathlib.Path
        Directorio que contiene los archivos FASTQ comprimidos (.fastq.gz / .fq.gz)
        de las lecturas a alinear.
    out_dir : str o pathlib.Path
        Directorio donde se escribirán los archivos BAM ordenados e indexados.
        Si no existe, se creará.
    threads : int, opcional
        Número de hilos que se usarán tanto en BWA como en samtools sort.
        Por defecto 8.
    sort_mem : str, opcional
        Cantidad de memoria por hilo para `samtools sort`, en formato
        aceptado por samtools (por ejemplo, "2G", "500M"). Por defecto "2G".
    tmpdir : str, opcional
        Directorio temporal donde `samtools sort` crea archivos temporales.
        Por defecto "/dev/shm".

    Salida
    ------
    No retorna ningún valor. Genera como efecto colateral:
      - Archivos *.sorted.bam en `out_dir`.
      - Archivos *.sorted.bam.bai (índices) en `out_dir`.

    Raises
    ------
    FileNotFoundError
        Si no se encuentra el archivo del genoma o no hay FASTQ en reads_dir.
    RuntimeError
        Si no se pueden identificar pares válidos R1/R2.
    
    """

    reads_dir = Path(reads_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genome_path = Path(bwa_index)
    if not genome_path.exists():
        raise FileNotFoundError(f"No existe el archivo del genoma: {genome_path}")

    # --- Descomprimir temporalmente si es .gz ---
    created_temp_fasta = False
    if genome_path.suffix == ".gz":
        decompressed_path = genome_path.with_suffix("")  # quita solo .gz
        if not decompressed_path.exists():
            print(f"Descomprimiendo genoma temporalmente: {genome_path.name}")
            with gzip.open(genome_path, "rb") as f_in, open(decompressed_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            created_temp_fasta = True
        genome_path = decompressed_path

    # --- Verificar/generar índice BWA ---
    index_files = [f"{genome_path}.{ext}" for ext in ["bwt", "pac", "ann", "amb", "sa"]]
    if any(not Path(f).exists() for f in index_files):
        print(f"Generando índice BWA para {genome_path.name}...")
        subprocess.run(["bwa", "index", str(genome_path)], check=True)
        print("Índice BWA creado correctamente.")

    # --- Buscar archivos FASTQ ---
    fastqs = sorted(list(reads_dir.glob("*.fastq.gz")) + list(reads_dir.glob("*.fq.gz")))
    if not fastqs:
        raise FileNotFoundError(f"No se encontraron archivos FASTQ en {reads_dir}")

    # --- Agrupar R1/R2 ---
    def sample_key(p: Path):
        
        name = p.name
        pair = None
        if re.search(r"(?:^|[_-])R?1(?:[_-]val[_-]?1)?(?=\.(?:fastq|fq)\.gz$)", name, re.IGNORECASE):
            pair = 1
        elif re.search(r"(?:^|[_-])R?2(?:[_-]val[_-]?2)?(?=\.(?:fastq|fq)\.gz$)", name, re.IGNORECASE):
            pair = 2
        base = re.sub(r"(?:^|[_-])R?[12](?:[_-]val[_-]?[12])?(?=\.(?:fastq|fq)\.gz$)", "", name, flags=re.IGNORECASE)
        return (base, pair)

    groups = {}
    for fq in fastqs:
        base, pair = sample_key(fq)
        if pair is None:
            print(f"[WARN] No se reconoce si es R1/R2: {fq.name} — se omite.")
            continue
        groups.setdefault(base, {})[pair] = fq

    pairs = [(d[1], d[2], base) for base, d in groups.items() if 1 in d and 2 in d]
    if not pairs:
        raise RuntimeError("No se encontraron pares R1/R2 válidos.")

    # --- Ejecutar alineamiento ---
    for r1, r2, base in pairs:
        sample_tag = Path(base).stem
        bam_out = out_dir / f"{sample_tag}.sorted.bam"
        if bam_out.exists():
            print(f"[INFO] {bam_out.name} ya existe — se omite.")
            continue

        print(f"\n Alineando muestra: {sample_tag}")
        print(f"   R1: {r1.name}")
        print(f"   R2: {r2.name}")
        print(f"   Hilos: {threads}, Memoria sort: {sort_mem}, TMP: {tmpdir}")

        p1 = subprocess.Popen(
            ["bwa", "mem", "-t", str(threads), str(genome_path), str(r1), str(r2)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p2 = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(
            [
                "samtools", "sort",    #Ordenar con samtools
                "-@", str(threads),
                "-m", sort_mem,
                "-T", f"{tmpdir}/{sample_tag}",
                "-o", str(bam_out),
                "-"
            ],
            stdin=p2.stdout
        )

        p1.stdout.close()
        p2.stdout.close()
        p3.communicate()

        subprocess.run(["samtools", "index", "-@", str(threads), str(bam_out)], check=True)
        print(f"{bam_out.name} y su índice generados correctamente.")

    # --- Eliminar genoma temporal ---
    if created_temp_fasta:
        print(f"Eliminando genoma temporal: {genome_path.name}")
        genome_path.unlink()

    print("\n Alineamientos completados correctamente.")

def generar_conteos_coverageBed(bam_dir, gff_file, out_dir=None):
    """
    Genera archivos .count.txt por muestra usando bedtools coverage (coverageBed).

    Flujo general
    -------------
    1. Recorrer todos los archivos BAM en `bam_dir`.
    2. Para cada BAM, ejecutar `bedtools coverage` con:
         -a <archivo GFF/BED> (regiones anotadas)
         -b <archivo BAM>    (alineamientos)
    3. Guardar la salida de bedtools coverage en un archivo:
         <sample_name>.count.txt en el directorio `out_dir`
       (o en `bam_dir` si out_dir no se especifica).

    Parámetros
    ----------
    bam_dir : str o pathlib.Path
        Directorio con los archivos .bam alineados.
    gff_file : str o pathlib.Path
        Archivo GFF/BED con las regiones anotadas (genes, CDS, exones, etc.).
    out_dir : str o pathlib.Path, opcional
        Directorio de salida para los archivos .count.txt.
        Si es None, se utilizan los mismos archivos .bam_dir como salida.

    Salida
    ------
    No retorna nada. Como efecto colateral, crea uno o varios archivos
    *.count.txt en out_dir.

    Raises
    ------
    FileNotFoundError
        Si no se encuentra el GFF/BED o si no hay BAMs en bam_dir.
    """
    import subprocess
    from pathlib import Path

    bam_dir = Path(bam_dir)
    out_dir = Path(out_dir) if out_dir else bam_dir   #Si no hay out_dir, manda los out_file al mismo count dir (bam_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gff_file = Path(gff_file)
    if not gff_file.exists():
        raise FileNotFoundError(f"No se encontró el archivo GFF: {gff_file}")

    bams = sorted(bam_dir.glob("*.bam"))
    if not bams:
        raise FileNotFoundError(f"No se encontraron archivos BAM en {bam_dir}")

    print(f"\n=== INICIANDO CUANTIFICACIÓN CON bedtools coverage ===")
    print(f"Anotación: {gff_file}")
    print(f"Total de muestras: {len(bams)}\n")

    for bam in bams:
        sample_name = bam.stem.replace(".fq", "").replace(".sorted", "")
        out_file = out_dir / f"{sample_name}.count.txt"

        if out_file.exists():
            print(f"[INFO] {out_file.name} ya existe — se omite.")
            continue

        print(f"Contando lecturas en: {bam.name}")
        cmd = [
            "bedtools", "coverage",
            "-a", str(gff_file),
            "-b", str(bam)
        ]

        with open(out_file, "w") as fout:
            result = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(f"[ERROR] Falló el conteo de {bam.name}:\n{result.stderr}")
        else:
            print(f"[OK] Conteo generado: {out_file.name}")

    print("\n Conteo de lecturas completado con éxito.\n")


def _etiqueta_muestra(filename: str):
    """
 
    Deduce el nombre de la muestra a partir del nombre del archivo .count.txt.

    Soporta patrones como:
      - GSM4099077.count.txt
      - SRR10192868.clean.sort.count.txt
      - sampleGSMxxxx.clean.sort.count.txt

    La idea es obtener una etiqueta de muestra limpia que luego se usará
    como nombre de columna en la tabla de cobertura.

    Parámetros
    ----------
    filename : str
        Nombre del archivo de conteo (sin ruta completa).

    Returns
    -------
    str
        Nombre de la muestra deducido a partir del nombre de archivo.
    
    """
    base = re.split(r"[._]", filename)[0]
    return re.sub(r"^sample", "", base, flags=re.IGNORECASE)

def _lee_y_normaliza_count(path: str):
    """
    
    Lee un archivo .count.txt generado por bedtools coverage y devuelve
    un DataFrame en formato "largo" (long format) con columnas:

        ['gene_id', 'gene_name', 'sample', 'count']

    Flujo general
    -------------
    1. Leer el archivo de conteo como tabla tabulada, tomando solo dos columnas:
       - Columna 8: información de anotación (gene_information).
       - Columna 9: conteo (coverage / número de lecturas).
    2. Extraer `gene_id` y `gene_name` de la cadena gene_information
       usando expresiones regulares sobre campos tipo:
         ID=<id>;Name=<name>;...
    3. Convertir la columna count a numérica y reemplazar valores no
       numéricos por 0.
    4. Agrupar por (gene_id, gene_name) y sumar los conteos, por si
       hay líneas duplicadas del mismo gen.
    5. Añadir la columna `sample` con la etiqueta de la muestra deducida
       a partir del nombre del archivo.

    Parámetros
    ----------
    path : str
        Ruta al archivo .count.txt a leer.

    Returns
    -------
    pandas.DataFrame
        DataFrame con columnas ['gene_id', 'gene_name', 'sample', 'count'].
    
    """
    p = Path(path)
    sample = _etiqueta_muestra(p.name)

    # Leer columnas 9 y 10
    df = pd.read_csv(
        p, sep="\t", header=None, usecols=[8, 9],
        names=["gene_information", "count"],
        dtype={"gene_information": "string"}
    )

    # Extraer ID y Name
    df["gene_id"] = df["gene_information"].str.extract(r"locus_tag=([^;]+)", expand=False)
    df["gene_name"] = df["gene_information"].str.extract(r"Name=([^;]+)", expand=False)
    # Si no hay Name, usar gene_id como gene_name
    df["gene_name"].fillna(df["gene_id"], inplace=True)
    # Asegurar que los conteos sean numéricos; reemplazar no numéricos por 0
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0)

    # Sumar por gene_id (por si hay duplicados)
    df = df.groupby("gene_id", as_index=False).agg(
    gene_name=("gene_name", "first"),
    count=("count", "sum"),
    )

    df["sample"] = sample       

    return df

def combinar_cobertura(path_in, out_csv, workers=None, pattern=None):
    """
    Combina archivos de conteo (*.count.txt) en una tabla de cobertura
    genes x muestras, lista para ser usada en DESeq2 u otros análisis.

    Flujo general
    -------------
    1. Localiza los archivos de conteo en `path_in`.
       - Por defecto, busca:
           *.count.txt
           *.clean.sort.count.txt
    2. Usa un `ProcessPoolExecutor` para leer en paralelo cada archivo
       con `_lee_y_normaliza_count()`, obteniendo tablas en formato largo.
    3. Concatena todos los DataFrames "largos" en uno solo.
    4. Realiza un pivot (pivot_table) para pasar a formato ancho:
       - Índice: (gene_id, gene_name)
       - Columnas: sample
       - Valores: count
    5. Ordena las columnas por nombre de muestra y guarda el resultado
       en `out_csv`.

    Parámetros
    ----------
    path_in : str o pathlib.Path
        Carpeta con los archivos de conteo (*.count.txt).
    out_csv : str o pathlib.Path
        Ruta del archivo CSV de salida donde se escribirá la tabla final
        de cobertura (genes x muestras).
    workers : int o None, opcional
        Número de procesos a utilizar para la lectura en paralelo.
        Si es None, se toma min(8, número de CPUs disponibles).
    pattern : str o None, opcional
        Patrón de búsqueda personalizado para los archivos de conteo.
        Si es None, se usan por defecto "*.count.txt" y "*.clean.sort.count.txt".

    Salida
    ------
    No retorna nada. Escribe en disco un archivo CSV con la tabla de
    cobertura en formato ancho.

    Raises
    ------
    FileNotFoundError
        Si no se encuentran archivos de conteo en path_in.
    RuntimeError
        Si no se pudo leer ningún archivo de conteo válido.
    """    
    path_in = Path(path_in)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # Selección de archivos de conteo según el patrón
    if pattern is None:
        files = sorted(path_in.glob("*.count.txt")) + sorted(path_in.glob("*.clean.sort.count.txt"))
    else:
        files = sorted(path_in.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No se encontraron archivos de conteo en {path_in}")
    
    # Determinar número de procesos a usar
    if workers is None:
        workers = min(8, max(1, os.cpu_count() or 2))

    print(f"[INFO] Combinando {len(files)} archivos de conteo usando {workers} procesos...")

    frames = []
    # Leer archivos de conteo en paralelo
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_lee_y_normaliza_count, str(f)): f for f in files}
        for fut in as_completed(futs):
            try:
                frames.append(fut.result())
            except Exception as e:
                print(f"[ERROR] Falló al leer {futs[fut]}: {e}")

    if not frames:
        raise RuntimeError("No se pudo leer ningún archivo de conteo válido.")
    
    # Concatenar en un DataFrame "largo"
    largo = pd.concat(frames, ignore_index=True)

    duplicados = largo.duplicated(subset=["gene_id", "sample"], keep=False)
    n_dup_pairs = duplicados.sum()
    if n_dup_pairs > 0:
        print(
            f"[WARN] Se encontraron {n_dup_pairs} filas duplicadas en (gene_id, sample) "
            "antes de colapsar. Se resolverán eligiendo un solo gene_name por gen."
        )

    largo["es_npyp"] = largo["gene_name"].str.startswith(("NP_", "YP_"), na=False)

    # Ordenar para que, para cada (gene_id, sample), queden primero los gene_name bonitos
    largo = largo.sort_values(["gene_id", "sample", "es_npyp"])

    # Para cada (gene_id, sample), quedarnos con la primera fila
    largo = largo.drop_duplicates(subset=["gene_id", "sample"], keep="first")

    # Quitar columna auxiliar
    largo = largo.drop(columns=["es_npyp"])

    #buscar duplicados
    if largo.duplicated(subset=["gene_id", "sample"]).any():
        raise RuntimeError(
            "Persisten duplicados en (gene_id, sample) después de colapsar. "
            "Revisar la lógica de combinación de cobertura."
        )

    # Convertir a formato ancho (genes en filas, muestras en columnas)
    ancho = largo.pivot_table(
        index=["gene_id", "gene_name"],
        columns="sample",
        values="count",
        aggfunc="sum",
        fill_value=0
    )

    ancho = ancho.reindex(sorted(ancho.columns), axis=1)
    ancho.reset_index().to_csv(out_csv, index=False)

    n_genes = ancho.index.get_level_values("gene_id").nunique()
    print(f"[OK] Tabla de cobertura generada correctamente: {out_csv}")
    print(f"[INFO] Genes únicos en tabla de cobertura: {n_genes}")

