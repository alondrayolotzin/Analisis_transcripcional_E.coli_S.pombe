"""
main_deseq_clasificar.py
Autor: Alondra Márquez
Fecha: 2025-11-24

Descripción
-----------
Script para:

1. Leer una tabla de resultados crudos de DESeq2 (salida de main_deseq_fit.py),
   que contiene columnas como baseMean, log2FoldChange, pvalue, padj, etc.
2. Clasificar los genes en tres categorías (UP, DOWN, Non-DE) de acuerdo con
   umbrales de padj y |log2FoldChange|.
3. Guardar una tabla de resultados con la nueva columna de clasificación
   (DE_status) y generar, opcionalmente, listas de genes UP y DOWN.

Este script sólo hace post-procesamiento a partir de la tabla ya generada.

Ejemplo de uso 
--------------
python main_deseq_clasificar.py \
    --resultados results/DESeq2_Ecoli/deseq2_resultados_crudos.csv \
    --sep "," \
    --lfc-col log2FoldChange \
    --padj-col padj \
    --padj-umbral 0.05 \
    --lfc-umbral 1.0 \
    --outdir results/DESeq2_Ecoli_clasificado \
    --log-file logs/deseq_clasificar_ecoli.log
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from modules.expresion_diff import (
    clasificar_genes,
    guardar_resultados,
)


def main():
    """
    Clasifica los resultados de DESeq2 en categorías UP/DOWN/Non-DE según
    umbrales de significancia (padj) y magnitud de cambio (|log2FoldChange|).

    Parámetros de entrada
    ---------------------
    La lee los argumentos de línea de comandos.

    --resultados : str
        Ruta al archivo CSV/TSV con los resultados crudos de DESeq2.
        Normalmente es la salida de `main_deseq_fit.py`, que contiene columnas
        como 'baseMean', 'log2FoldChange', 'pvalue', 'padj', etc.
    --sep : str, opcional
        Separador del archivo de entrada de resultados (por defecto ",").
    --lfc-col : str, opcional
        Nombre de la columna que contiene el log2FoldChange. Por defecto
        'log2FoldChange'.
    --padj-col : str, opcional
        Nombre de la columna que contiene el p-valor ajustado (padj).
        Por defecto 'padj'.
    --padj-umbral : float, opcional
        Umbral máximo de padj para considerar un gen como diferencialmente
        expresado (DE). Por defecto 0.05.
    --lfc-umbral : float, opcional
        Umbral mínimo de la magnitud del log2FoldChange (|log2FC|) para
        considerar un gen como DE. Por defecto 1.0.
    --outdir : str
        Directorio de salida donde se guardarán:
        - la tabla de resultados clasificados (con columna 'DE_status'),
        - y, opcionalmente, las listas de genes UP y DOWN.
    --log-file : str, opcional
        Ruta al archivo de log. Si no se especifica, se usa:
        <outdir>/deseq_clasificar.log

    Salida
    ------

    * Crea el directorio de salida especificado en --outdir (si no existe).
    * Escribe un archivo de log con mensajes de progreso y posibles errores.
    * Genera los siguientes archivos en el directorio de salida:

        - deseq2_resultados_clasificados.csv
          Tabla de resultados de DESeq2 con una nueva columna 'DE_status' que
          clasifica cada gen como 'UP', 'DOWN' o 'Non-DE'.

        - genes_UP.txt   (opcional)
          Lista de IDs de genes clasificados como 'UP', una por línea.

        - genes_DOWN.txt (opcional)
          Lista de IDs de genes clasificados como 'DOWN', una por línea.

      Las listas UP/DOWN sólo se generan si la tabla clasificada contiene
      una columna llamada 'gene_id'. En caso contrario, se emite una advertencia
      en el log.

    """
    # === Definición y parseo de argumentos de línea de comandos ===
    parser = argparse.ArgumentParser(
        description=(
            "Clasificación de resultados de DESeq2 en UP/DOWN/Non-DE según umbrales.\n\n"
            "Lee una tabla de resultados crudos (salida de main_deseq_fit.py) y aplica\n"
            "umbrales de padj y |log2FoldChange| para etiquetar cada gen."
        )
    )
    parser.add_argument("--resultados",required=True,help="CSV/TSV con resultados crudos de DESeq2 (salida de main_deseq_fit.py).",)
    parser.add_argument("--sep",default=",",help="Separador del archivo de resultados (por defecto ',').",)
    parser.add_argument("--lfc-col", default="log2FoldChange",help="Nombre de la columna con log2FC (por defecto 'log2FoldChange').",)
    parser.add_argument("--padj-col",default="padj", help="Nombre de la columna con padj (por defecto 'padj').",)
    parser.add_argument("--padj-umbral",type=float,default=0.05,help="Umbral máximo de padj para considerar un gen DE (por defecto 0.05).",)
    parser.add_argument("--lfc-umbral",type=float,default=1.0,help="Umbral mínimo de |log2FC| para considerar un gen DE (por defecto 1.0).",)
    parser.add_argument("--outdir",required=True,help="Directorio de salida para resultados clasificados y archivo de log.",)
    parser.add_argument("--log-file",default=None,
        help=(
            "Ruta al archivo de log (por defecto: <outdir>/deseq_clasificar.log). "
            "El directorio se creará si no existe."
        ),
    )

    args = parser.parse_args()

    # === Preparación de directorio de salida ===
    directorio_salida = Path(args.outdir).resolve()
    directorio_salida.mkdir(parents=True, exist_ok=True)

    # === Configuración de logging (estilo main_download.py) ===
    ruta_log = Path(args.log_file) if args.log_file else directorio_salida / "deseq_clasificar.log"
    log_dir = ruta_log.parent
    log_dir.mkdir(parents=True, exist_ok=True)

    log_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    manejador_archivo = logging.FileHandler(ruta_log, mode="w")
    manejador_archivo.setFormatter(log_formatter)

    manejador_consola = logging.StreamHandler(sys.stdout)
    manejador_consola.setFormatter(log_formatter)

    logging.getLogger().handlers = []  # limpia handlers previos
    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(manejador_archivo)
    logging.getLogger().addHandler(manejador_consola)

    logger = logging.getLogger(__name__)

    logger.info("=== Inicio de la clasificación de genes DE (UP/DOWN/Non-DE) ===")

    # 1) Leer resultados crudos de DESeq2
    ruta_resultados = Path(args.resultados)
    if not ruta_resultados.exists():
        logger.error(f"No se encontró el archivo de resultados: {ruta_resultados}")
        sys.exit(1)

    tabla_resultados_crudos = pd.read_csv(ruta_resultados, sep=args.sep)
    logger.info(
        f"Resultados crudos leídos desde {ruta_resultados} "
        f"({tabla_resultados_crudos.shape[0]} filas, "
        f"{tabla_resultados_crudos.shape[1]} columnas)."
    )

    # 2) Clasificar genes según padj y log2FC
    tabla_clasificada, resumen_categorias = clasificar_genes(
        res_df=tabla_resultados_crudos,
        lfc_col=args.lfc_col,
        padj_col=args.padj_col,
        padj_umbral=args.padj_umbral,
        lfc_umbral=args.lfc_umbral,
        nueva_columna="DE_status",  # puedes cambiar a "Expresion" si quieres
    )

    # 3) Guardar tabla de resultados clasificados
    ruta_resultados_clasificados = directorio_salida / "deseq2_resultados_clasificados.csv"
    guardar_resultados(
        df=tabla_clasificada,
        path_salida=ruta_resultados_clasificados,
        sep=",",
        index=True,
    )

    # 4) (Opcional) guardar listas UP y DOWN separadas, si hay gene_id
    ruta_genes_up = directorio_salida / "genes_UP.txt"
    ruta_genes_down = directorio_salida / "genes_DOWN.txt"

    if "DE_status" in tabla_clasificada.columns and "gene_id" in tabla_clasificada.columns:
        tabla_clasificada.loc[
            tabla_clasificada["DE_status"] == "UP", "gene_id"
        ].to_csv(ruta_genes_up, index=False, header=False)

        tabla_clasificada.loc[
            tabla_clasificada["DE_status"] == "DOWN", "gene_id"
        ].to_csv(ruta_genes_down, index=False, header=False)

        logger.info(f"Lista de genes UP guardada en:   {ruta_genes_up}")
        logger.info(f"Lista de genes DOWN guardada en: {ruta_genes_down}")
    else:
        logger.warning(
            "No se encontraron columnas 'DE_status' y/o 'gene_id'; "
            "no se generaron listas de genes UP/DOWN."
        )

    logger.info(f"Resumen de categorías DE: {resumen_categorias}")
    logger.info(f"Resultados clasificados guardados en: {ruta_resultados_clasificados}")
    logger.info("=== Clasificación finalizada correctamente ===")
    print("\nClasificación de genes DE completada. Revisar logs y resultados para más detalles.\n")


if __name__ == "__main__":
    main()
