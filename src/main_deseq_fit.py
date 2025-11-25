"""
main_deseq_fit.py
Autor: Alondra Márquez
Fecha: 2025-11-24

Descripción
-----------
Pipeline para:

1. Leer la tabla de cobertura generada a partir de coverageBed/combinar_cobertura
   (formato típico: gene_id, gene_name, muestras...).
2. Limpiar la fila del cromosoma completo (ANONYMOUS), si está presente.
3. Construir la matriz de conteos (genes x muestras).
4. Renombrar las columnas de la matriz de conteos usando un archivo de mapeo
   (sample_id,nuevo_nombre), por ejemplo:
       GSM4099077 -> Gluc_eco_rep1
5. Usar un archivo de metadata (sample_id,condition,medium,organism,replicate)
   para definir la condición de cada muestra (Gluc vs Glyc).
6. Aplicar un filtro de expresión mínima basado en CPM.
7. Construir el objeto DeseqDataSet (PyDESeq2) usando la columna 'condition'
   de la metadata.
8. Ejecutar el análisis de expresión diferencial con DESeq2.
9. Guardar una tabla de resultados crudos (baseMean, log2FoldChange, pvalue, padj, etc.)
   lista para gráficas y pasos posteriores.

Ejemplo de uso 
--------------
python main_deseq_fit.py \
    --counts ../results/tabla_cobertura_Ecoli.csv \
    --sep "," \
    --mapeo ../data/E.coli_samples.csv \
    --metadata ../data/metadata_ecoli_defseq.csv \
    --condicion-prueba Gluc \
    --condicion-control Glyc \
    --umbral-cpm 5 \
    --min-muestras 3 \
    --outdir results/DE/DESeq2_Ecoli \
    --log logs/deseq_fit_ecoli.log
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

from modules.expresion_diff import (
    cambiar_nombre_col,
    filtrar_por_cpm,
    crear_objeto_dds,
    deseq2,
    estadisticas,
    guardar_resultados,
)


def main():
    """
    Ejecuta parte del pipeline de ajuste del modelo de expresión diferencial
    usando PyDESeq2 a partir de una tabla de cobertura.

    Parámetros de entrada
    ---------------------
    La función recibe parámetros desde la línea de comandos, los principales son: 

    --counts : str
        Ruta al archivo CSV/TSV con la tabla de cobertura. Debe incluir las
        columnas 'gene_id', 'gene_name' y una columna por muestra con conteos.
    --sep : str, opcional
        Separador del archivo de conteos (por defecto ",").
    --mapeo : str
        Ruta a un CSV con el mapeo de nombres de muestras. Debe contener
        al menos las columnas 'sample_id' y 'nuevo_nombre'.
    --metadata : str
        Ruta al CSV de metadata con información de cada muestra. Debe
        contener al menos las columnas:
        'sample_id', 'condition', 'medium', 'organism', 'replicate'.
    --umbral-cpm : float, opcional
        Umbral mínimo de CPM para considerar que un gen está expresado en
        una muestra (por defecto 1.0).
    --min-muestras : int, opcional
        Número mínimo de muestras que deben superar el umbral de CPM para
        conservar un gen (por defecto 2).
    --condicion-prueba : str
        Nombre de la condición de prueba (tratamiento) en la columna
        'condition' de la metadata (por ejemplo 'Glyc').
    --condicion-control : str
        Nombre de la condición de referencia (control) en la columna
        'condition' de la metadata (por ejemplo 'Gluc').
    --outdir : str
        Directorio de salida donde se guardarán los archivos resultantes.
    --log : str, opcional
        Ruta al archivo de log. Si no se especifica, se usa
        <outdir>/deseq_fit.log.

    Salida
    ------

    * Crea el directorio de salida especificado en --outdir (si no existe).
    * Escribe un archivo de log con mensajes de progreso y posibles errores.
    * Genera los siguientes archivos en el directorio de salida:

        - counts_filtrados_cpm.csv
          Matriz de conteos filtrada por CPM (genes x muestras).

        - deseq2_resultados_crudos.csv
          Tabla de resultados crudos del análisis de expresión diferencial
          (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, etc.).

    Además, muestra en pantalla un mensaje final indicando que el ajuste
    DESeq2 se completó correctamente, o termina el programa con código de
    salida 1 si ocurre algún error crítico (por ejemplo, archivos faltantes
    o formatos incompatibles).
    """

    parser = argparse.ArgumentParser(description=(
            "Ajuste de modelo de expresión diferencial con PyDESeq2.\n"
            "Flujo: tabla de cobertura -> limpiar / renombrar -> "
            "usar metadata -> filtro CPM -> DESeq2 -> guardar resultados crudos."
        )
    )
    parser.add_argument("--counts",required=True, help="Ruta al CSV/TSV con la tabla de cobertura (gene_id,gene_name,muestras).",)
    parser.add_argument("--sep",default=",",help="Separador del archivo de conteos (por defecto ',').",)
    parser.add_argument("--mapeo",required=True, help=("CSV con mapeo de nombres de muestras. ""Debe contener las columnas 'sample_id' y 'nuevo_nombre'."),)
    parser.add_argument("--metadata",required=True,help=("CSV de metadata con columnas 'sample_id' y 'condition' " "(además de medium, organism, replicate)."))

    parser.add_argument("--umbral-cpm",type=float,default=1.0,help="Umbral mínimo de CPM (por defecto 1.0).",)
    parser.add_argument("--min-muestras",type=int,default=2,help="Número mínimo de muestras con CPM >= umbral (por defecto 2).",)

    parser.add_argument("--condicion-prueba",required=True,help="Nombre de la condición de prueba (tratamiento), por ejemplo 'Glyc'.",)
    parser.add_argument("--condicion-control",required=True,help="Nombre de la condición de referencia (control), por ejemplo 'Gluc'.",)

    # Salida y log
    parser.add_argument("--outdir",required=True,help="Directorio de salida para resultados.",)
    parser.add_argument("--log",default=None,help="Ruta al archivo de log (por defecto: <outdir>/deseq_fit.log).",
)

    args = parser.parse_args()

    # === Preparación de directorios de salida ===
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    log_path = Path(args.log) if args.log else outdir / "deseq_fit.log"
    log_dir = os.path.dirname(log_path)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    # === Configuración de logging ===
    log_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(log_formatter)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)

    logging.getLogger().handlers = []  # limpia cualquier handler previo
    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(file_handler)
    logging.getLogger().addHandler(console_handler)

    logger = logging.getLogger(__name__)

    logger.info("=== Iniciando ajuste DESeq2 (PyDESeq2) ===")

    # Leer tabla de cobertura
    counts_path = Path(args.counts)
    if not counts_path.exists():
        logger.error(f"No se encontró el archivo de conteos: {counts_path}")
        sys.exit(1)

    df_cobertura = pd.read_csv(counts_path, sep=args.sep)
    logger.info(
        f"Tabla de cobertura leída desde {counts_path} "
        f"({df_cobertura.shape[0]} filas, {df_cobertura.shape[1]} columnas)."
    )

    # 2) Limpiar fila del cromosoma completo (ANONYMOUS), si existe
    if {"gene_id", "gene_name"}.issubset(df_cobertura.columns):
        mascara_cromosoma = (
            df_cobertura["gene_id"].astype(str).str.contains(":", regex=False)
            & (df_cobertura["gene_name"] == "ANONYMOUS")
        )
        filas_antes = df_cobertura.shape[0]
        df_cobertura = df_cobertura.loc[~mascara_cromosoma].copy()
        filas_despues = df_cobertura.shape[0]
        if filas_despues < filas_antes:
            logger.info(
                f"Se eliminaron {filas_antes - filas_despues} filas correspondientes "
                "al cromosoma completo (ANONYMOUS)."
            )

        # Construir matriz de conteos: índice = gene_id, columnas = muestras
        df_counts = df_cobertura.set_index("gene_id").drop(columns=["gene_name"])
    else:
        logger.error(
            "La tabla de cobertura no tiene columnas 'gene_id' y 'gene_name'. "
            "Revisa el formato de entrada."
        )
        sys.exit(1)

    n_genes = df_counts.shape[0]
    n_unique = df_counts.index.nunique()
    logger.info(
        f"Matriz de conteos inicial: {n_genes} filas; gene_id únicos: {n_unique}"
    )

    if n_genes != n_unique:
        logger.warning(
        "Se detectaron gene_id duplicados en la matriz de conteos. "
        "Revisa la tabla de cobertura o el módulo de cobertura."
    )

    logger.info(
        f"Matriz de conteos inicial: {df_counts.shape[0]} genes x "
        f"{df_counts.shape[1]} muestras."
    )

    # 3) Renombrar columnas con archivo de mapeo (sample_id a nuevo_nombre)
    mapeo_path = Path(args.mapeo)
    if not mapeo_path.exists():
        logger.error(f"No se encontró el archivo de mapeo: {mapeo_path}")
        sys.exit(1)

    df_counts = cambiar_nombre_col(
        df_counts,
        mapeo_path,
        "sample_id",
        "nuevo_nombre",
    )
    logger.info(
        f"Columnas renombradas. Muestras actuales: {list(df_counts.columns)}"
    )

    # 4) Leer metadata y alinearla con la matriz de conteos
    metadata_path = Path(args.metadata)
    if not metadata_path.exists():
        logger.error(f"No se encontró el archivo de metadata: {metadata_path}")
        sys.exit(1)

    metadata = pd.read_csv(metadata_path)
    columnas_esperadas = {"sample_id", "condition", "medium", "organism", "replicate"}
    if not columnas_esperadas.issubset(metadata.columns):
        logger.error(
            "La metadata no tiene las columnas esperadas "
            f"{columnas_esperadas}. Columnas encontradas: {list(metadata.columns)}"
        )
        sys.exit(1)

    metadata = metadata.set_index("sample_id")

    # Asegurar que metadata y conteos se refieren a las mismas muestras
    muestras_conteos = set(df_counts.columns)
    muestras_metadata = set(metadata.index)

    muestras_comunes = sorted(muestras_conteos & muestras_metadata)
    if len(muestras_comunes) == 0:
        logger.error(
            "No hay intersección entre las muestras de la matriz de conteos "
            "y las de la metadata."
        )
        sys.exit(1)

    # Reordenar y filtrar para quedarse sólo con las muestras comunes
    df_counts = df_counts.loc[:, muestras_comunes]
    metadata = metadata.loc[muestras_comunes]

    logger.info(
        f"Tras alinear con metadata: {df_counts.shape[0]} genes x "
        f"{df_counts.shape[1]} muestras."
    )
    logger.info(
        f"Condiciones presentes en metadata: {metadata['condition'].unique()}"
    )

    # 5) Filtro por CPM
    path_counts_filtrados = outdir / "counts_filtrados_cpm.csv"
    df_filtrada = filtrar_por_cpm(
        df_counts,
        args.umbral_cpm,
        args.min_muestras,
        path_counts_filtrados,
    )

    logger.info(
        f"Matriz filtrada por CPM: {df_filtrada.shape[0]} genes x "
        f"{df_filtrada.shape[1]} muestras."
    )

    # 6) Crear DeseqDataSet a partir de conteos transpuestos y metadata
    matriz_conteos_transpuesta = df_filtrada.T
    dds = crear_objeto_dds(
        matriz_conteos_transpuesta,
        metadata,
        "condition",
    )
    logger.info("Objeto DeseqDataSet creado correctamente.")

    # 7) Ejecutar DESeq2
    dds = deseq2(dds)

    # 8) Estadísticas de expresión diferencial
    res_df = estadisticas(dds,args.condicion_prueba,args.condicion_control,"condition",)

    # 9) Guardar resultados crudos
    path_res_raw = outdir / "deseq2_resultados_crudos.csv"
    guardar_resultados(res_df,path_res_raw,",",True,)

    logger.info(f"Resultados crudos de DESeq2 guardados en: {path_res_raw}")
    logger.info("=== Ajuste DESeq2 finalizado correctamente ===")
    print("\nAjuste DESeq2 completado. Revisar logs y resultados para más detalles.\n")


if __name__ == "__main__":
    main()
