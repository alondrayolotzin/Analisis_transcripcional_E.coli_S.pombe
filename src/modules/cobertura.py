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
    Alinea lecturas pareadas (R1/R2) con BWA MEM y produce archivos BAM ordenados e indexados.
    - Acepta referencia .fa/.fna o comprimida (.fa.gz, .fna.gz)
    - Descomprime temporalmente si es necesario
    - Usa samtools con parámetros configurables
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

        print(f"\n🧬 Alineando muestra: {sample_tag}")
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

    Parámetros
    ----------
    bam_dir : str | Path
        Directorio con los archivos .bam alineados.
    gff_file : str | Path
        Archivo GFF/BED con las regiones anotadas.
    out_dir : str | Path, opcional
        Directorio de salida para los .count.txt (por defecto, usa el mismo bam_dir).
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


def _etiqueta_muestra(filename: str) -> str:
    """
    Deduce el nombre de la muestra a partir del nombre del archivo .count.txt.
    Soporta:
      - GSM4099077.count.txt
      - SRR10192868.clean.sort.count.txt
      - sampleGSMxxxx.clean.sort.count.txt
    """
    base = re.split(r"[._]", filename)[0]
    return re.sub(r"^sample", "", base, flags=re.IGNORECASE)

def _lee_y_normaliza_count(path: str):
    """
    Lee un archivo .count.txt y devuelve un DataFrame con:
    ['gene_id', 'gene_name', 'sample', 'count']
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
    df["gene_id"] = df["gene_information"].str.extract(r"ID=([^;]+)", expand=False)
    df["gene_name"] = df["gene_information"].str.extract(r"Name=([^;]+)", expand=False)
    df["gene_name"].fillna(df["gene_id"], inplace=True)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0)

    # Sumar por gene_id (por si hay duplicados)
    df = df.groupby(["gene_id", "gene_name"], as_index=False)["count"].sum()
    df["sample"] = sample

    return df

def combinar_cobertura(path_in, out_csv, workers=None, pattern=None):
    """
    Combina archivos *.count.txt en una tabla de cobertura en formato ancho (genes x muestras).
    Usa procesamiento en paralelo para acelerar la lectura y combinación.

    Parámetros
    ----------
    path_in : str | Path
        Carpeta con los archivos de conteo.
    out_csv : str | Path
        Ruta del CSV de salida.
    workers : int | None
        Número de procesos (por defecto: min(8, núm. de CPUs)).
    pattern : str | None
        Patrón opcional (por defecto: *.count.txt y *.clean.sort.count.txt)
    """
    path_in = Path(path_in)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if pattern is None:
        files = sorted(path_in.glob("*.count.txt")) + sorted(path_in.glob("*.clean.sort.count.txt"))
    else:
        files = sorted(path_in.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No se encontraron archivos de conteo en {path_in}")

    if workers is None:
        workers = min(8, max(1, os.cpu_count() or 2))

    print(f"[INFO] Combinando {len(files)} archivos de conteo usando {workers} procesos...")

    frames = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_lee_y_normaliza_count, str(f)): f for f in files}
        for fut in as_completed(futs):
            try:
                frames.append(fut.result())
            except Exception as e:
                print(f"[ERROR] Falló al leer {futs[fut]}: {e}")

    if not frames:
        raise RuntimeError("No se pudo leer ningún archivo de conteo válido.")

    largo = pd.concat(frames, ignore_index=True)

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

    print(f"[OK] Tabla de cobertura generada correctamente: {out_csv}")

