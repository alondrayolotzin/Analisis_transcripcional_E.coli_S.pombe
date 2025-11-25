"""
main_aligments_counts.py
Autor: Miryam Zamora
Fecha: 2025-11-23

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
python src/main_aligments_counts.py --accion combinar \
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

    python scr/main_aligments_counts.py --accion todo \
  --genomas_dir data/Escherichia_coli_GCF_000005845.2 \
  --genoma GCF_000005845.2_ASM584v2_genomic.fna.gz \
  --gff data/Escherichia_coli_GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.gff.gz \
  --reads_dir results/Cleaned_fastq/Escherichia_coli_str._K-12_substr._MG1655 \
  --out_dir data/BWA/e_coli_bam \
  --count_dir data/BWA/e_coli_bam \
  --out_csv results/tabla_cobertura_Ecoli.csv \
  --threads 8
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
    """
    Función principal del script de cobertura.

    Se encarga de:
      1. Parsear los argumentos de línea de comandos.
      2. Según el valor de --accion, ejecutar una o varias fases del pipeline:
         - 'alinear'  : sólo alineamiento con BWA.
         - 'contar'   : sólo conteo de lecturas por gen.
         - 'combinar' : sólo combinación de archivos *.count.txt en una tabla final.
         - 'todo'     : ejecuta alinear + contar + combinar en secuencia.
      3. Verificar la existencia de archivos y directorios necesarios.
      4. Llamar a las funciones del módulo `modules.cobertura` que realizan el trabajo:
         - alinear_con_bwa()
         - generar_conteos_coverageBed()
         - combinar_cobertura()

    Parámetros
    ----------
    lee todos los parámetros desde la línea de comandos mediante argparse.

    Salida
    ------
    No retorna ningún valor. Si no hay errores, genera como efecto
    colateral:
      - Archivos BAM/SAM alineados (fase 'alinear').
      - Archivos de conteo por gen (*.count.txt) (fase 'contar').
      - Un archivo CSV final con la tabla de cobertura (fase 'combinar').
    Si ocurre un error crítico (por ejemplo, falta un archivo), el
    programa termina con sys.exit() mostrando un mensaje descriptivo.
    """


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

    # ==============================================================
    # ACCIÓN: ALINEAR o TODO
    # ==============================================================

    if args.accion in ("alinear", "todo"):
        # Verificación de argumentos requeridos para el alineamiento
        if not args.genomas_dir or not args.genoma or not args.reads_dir:
            sys.exit("Error: --genomas_dir, --genoma y --reads_dir son obligatorios para 'alinear' o 'todo'.")

    # Construcción de rutas absolutas a partir de la raíz del proyecto
        genomas_dir = (ROOT / args.genomas_dir).resolve()
        genoma_fasta = (genomas_dir / args.genoma).resolve()
        reads_dir = (ROOT / args.reads_dir).resolve()
        out_dir = (ROOT / args.out_dir).resolve()

    # Comprobación de existencia de genoma y directorio de lecturas
        if not genoma_fasta.exists():
            sys.exit(f"No se encontró el archivo del genoma: {genoma_fasta}")
        if not reads_dir.exists():
            sys.exit(f"No se encontró el directorio de lecturas: {reads_dir}")

        print("\n=== INICIANDO ALINEAMIENTO CON BWA ===\n")
        # Llama a la función de modules.cobertura que:
        #   - Indexa el genoma (si es necesario).
        #   - Alinea los FASTQ.
        #   - Ordena y genera archivos BAM finales.

        alinear_con_bwa(
            bwa_index=str(genoma_fasta),
            reads_dir=reads_dir,
            out_dir=out_dir,
            threads=args.threads,
            sort_mem=args.sort_mem,
            tmpdir=args.tmpdir,
        )

    # ==============================================================
    # ACCIÓN: CONTAR o TODO
    # ==============================================================

    if args.accion in ("contar", "todo"):
        # Verificación de argumentos requeridos para el conteo
        if not args.count_dir or not args.gff:
            sys.exit("Error: --count_dir y --gff son obligatorios para 'contar' o 'todo'.")

        count_dir = (ROOT / args.count_dir).resolve()
        gff_path = (ROOT / args.gff).resolve()

        # Comprobación de existencia de directorio de BAM y archivo GFF
        if not count_dir.exists():
            sys.exit(f"No se encontró el directorio de conteo: {count_dir}")
        if not gff_path.exists():
            sys.exit(f"No se encontró el archivo de anotaciones GFF/BED: {gff_path}")

        print("\n=== INICIANDO CONTEO DE LECTURAS CON bedtools coverage ===\n")
        # Llama a la función que recorre los BAM en count_dir y genera
        # archivos de conteo por gen/región anotada

        generar_conteos_coverageBed(count_dir, gff_path)

    # ==============================================================
    # ACCIÓN: COMBINAR o TODO
    # ==============================================================

    if args.accion in ("combinar", "todo"):
        # Verificación de argumentos requeridos para combinar
        if not args.count_dir or not args.out_csv:
            sys.exit("Error: --count_dir y --out_csv son obligatorios para 'combinar' o 'todo'.")

        count_dir = (ROOT / args.count_dir).resolve()
        out_csv = (ROOT / args.out_csv).resolve()
        out_csv.parent.mkdir(parents=True, exist_ok=True)

        print("\n=== GENERANDO TABLA DE COBERTURA FINAL ===\n")

        # Llama a la función que:
        #   - Lee todos los *.count.txt en count_dir.
        #   - Los combina en una sola tabla (genes x muestras).
        #   - Guarda el resultado en formato CSV

        combinar_cobertura(count_dir, out_csv, workers= args.threads)

    print("\n Proceso finalizado con éxito.\n")


if __name__ == "__main__":
    main()
