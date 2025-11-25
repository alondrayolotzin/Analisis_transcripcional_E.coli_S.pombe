#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main.cobertura.py
=================
Pipeline completo para:
  1) Alinear lecturas FASTQ a un genoma de referencia con BWA.
  2) Contar lecturas alineadas por gen usando featureCounts.
  3) Combinar archivos *.count.txt en una tabla final de cobertura.

Ejemplo de uso:
---------------
# Solo alinear
python src/main.cobertura.py --accion alinear \
    --genomas_dir genomas \
    --genoma GCF_000002945.2_ASM294v3_genomic.fna.gz \
    --reads_dir data/fastq \
    --out_dir data/BWA

# Solo contar
python src/main.cobertura.py --accion contar \
    --count_dir data/BWA \
    --gff genomas/GCF_000002945.2_ASM294v3_genomic.gff.gz \
    --threads 8

# Solo combinar
python src/main.cobertura.py --accion combinar \
    --count_dir data/BWA \
    --out_csv results/tabla_cobertura_Spombe.csv

# Pipeline completo (alinear + contar + combinar)
python src/main.cobertura.py --accion todo \
    --genomas_dir genomas \
    --genoma GCF_000002945.2_ASM294v3_genomic.fna.gz \
    --gff genomas/GCF_000002945.2_ASM294v3_genomic.gff.gz \
    --reads_dir data/fastq \
    --out_dir data/BWA \
    --count_dir data/BWA \
    --out_csv results/tabla_cobertura_Spombe.csv
"""
#nohup python src/main.cobertura.py --accion combinar --count_dir data/BWA/e_coli_bam/ --out_csv results/tabla_cobertura_Ecoli.csv --threads 8 > combinar_Ecoli.out 2>&1

import argparse
import sys
from pathlib import Path

# Ajustar la ruta para permitir importar desde src/modules
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))


from modules.cobertura import (
    alinear_con_bwa,
    generar_conteos_coverageBed,
    combinar_cobertura,
)


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline para alineamiento, conteo y combinación de archivos de cobertura."
    )

    parser.add_argument("--accion", choices=["alinear", "contar", "combinar", "todo"], required=True,
                        help="Acción a realizar: alinear | contar | combinar | todo")

    parser.add_argument("--genomas_dir", help="Directorio donde se encuentran los genomas de referencia (.fna.gz).")
    parser.add_argument("--genoma", help="Archivo del genoma de referencia (FASTA o .fna.gz).")
    parser.add_argument("--gff", help="Archivo GFF/GTF para featureCounts (requerido para 'contar' o 'todo').")
    parser.add_argument("--reads_dir", help="Directorio con archivos .fastq.gz (lecturas crudas o limpias).")
    parser.add_argument("--out_dir", default="data/BWA", help="Directorio de salida para alineamientos y conteos.")
    parser.add_argument("--count_dir", help="Directorio con archivos *.count.txt (para combinar o usar como salida de conteo).")
    parser.add_argument("--out_csv", help="Ruta de salida del archivo CSV con la tabla combinada.")

    parser.add_argument("--threads", type=int, default=8, help="Número de hilos a usar (por defecto: 8)")
    parser.add_argument("--sort-mem", type=str, default="2G", help="Memoria por hilo para Samtools sort (por defecto: 2G)")
    parser.add_argument("--tmpdir", default="/dev/shm", help="Directorio temporal para Samtools sort (por defecto: /dev/shm)")

    args = parser.parse_args()
    ROOT = Path(__file__).resolve().parent.parent

    # ========== ACCIÓN: ALINEAR o TODO ==========
    if args.accion in ("alinear", "todo"):
        if not args.genomas_dir or not args.genoma or not args.reads_dir:
            sys.exit("Error: --genomas_dir, --genoma y --reads_dir son obligatorios para 'alinear' o 'todo'.")

        genomas_dir = (ROOT / args.genomas_dir).resolve()
        genoma_fasta = (genomas_dir / args.genoma).resolve()
        reads_dir = (ROOT / args.reads_dir).resolve()
        out_dir = (ROOT / args.out_dir).resolve()

        if not genoma_fasta.exists():
            sys.exit(f"No se encontró el archivo del genoma: {genoma_fasta}")
        if not reads_dir.exists():
            sys.exit(f"No se encontró el directorio de lecturas: {reads_dir}")

        print("\n=== 🧬 INICIANDO ALINEAMIENTO CON BWA ===\n")
        alinear_con_bwa(
            bwa_index=str(genoma_fasta),
            reads_dir=reads_dir,
            out_dir=out_dir,
            threads=args.threads,
            sort_mem=args.sort_mem,
            tmpdir=args.tmpdir,
        )

    # ========== ACCIÓN: CONTAR o TODO ==========
    if args.accion in ("contar", "todo"):
        if not args.count_dir or not args.gff:
            sys.exit("Error: --count_dir y --gff son obligatorios para 'contar' o 'todo'.")

        count_dir = (ROOT / args.count_dir).resolve()
        gff_path = (ROOT / args.gff).resolve()

        if not count_dir.exists():
            sys.exit(f"No se encontró el directorio de conteo: {count_dir}")
        if not gff_path.exists():
            sys.exit(f"No se encontró el archivo de anotaciones GFF/BED: {gff_path}")

        print("\n=== INICIANDO CONTEO DE LECTURAS CON bedtools coverage ===\n")
        generar_conteos_coverageBed(count_dir, gff_path)

    # ========== ACCIÓN: COMBINAR o TODO ==========
    if args.accion in ("combinar", "todo"):
        if not args.count_dir or not args.out_csv:
            sys.exit("Error: --count_dir y --out_csv son obligatorios para 'combinar' o 'todo'.")

        count_dir = (ROOT / args.count_dir).resolve()
        out_csv = (ROOT / args.out_csv).resolve()
        out_csv.parent.mkdir(parents=True, exist_ok=True)

        print("\n=== GENERANDO TABLA DE COBERTURA FINAL ===\n")
        combinar_cobertura(count_dir, out_csv, workers= args.threads)

    print("\n Proceso finalizado con éxito.\n")


if __name__ == "__main__":
    main()
