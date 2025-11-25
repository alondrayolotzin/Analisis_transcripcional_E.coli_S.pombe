"""
expresion_diff.py
Autor: Alondra Márquez
Fecha: 2025-11-24

Descripción general
-------------------
Este módulo contiene las funciones centrales utilizadas en el proyecto para 
realizar análisis de expresión diferencial con PyDESeq2. Incluye rutinas para:

    • Renombrar columnas de la matriz de conteos usando un archivo de mapeo.
    • Calcular CPM (counts per million) y aplicar un filtro basado en expresión mínima.
    • Construir objetos DeseqDataSet a partir de conteos filtrados y metadata.
    • Ejecutar el flujo DESeq2 (estimación de size factors, dispersions y fitting del modelo).
    • Obtener estadísticas de expresión diferencial para un contraste específico.
    • Clasificar genes en categorías UP / DOWN / Non-DE según log2FC y padj.
    • Guardar tablas de resultados intermedios o finales.

Este módulo es utilizado directamente por los scripts:

    - main_deseq_fit.py
        Encargado de preparar los datos, crear el objeto DESeq2, 
        ejecutar el análisis y generar la tabla de resultados crudos.

    - main_deseq_clasificar.py
        Encargado de tomar la tabla generada por el script anterior y 
        clasificar los genes según criterios estadísticos establecidos.
"""
import logging
from pathlib import Path
from typing import Tuple, Dict
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

logger = logging.getLogger(__name__)



#  Cambiar nombres de columnas según archivo de mapeo


def cambiar_nombre_col(matriz_conteos: pd.DataFrame,ruta_mapeo: str | Path, columna_matriz: str,columna_nuevo: str,):
    """
    Cambia los nombres de las columnas de la matriz de conteos usando un CSV
    de mapeo, por ejemplo:

        sample_id,nuevo_nombre
        GSM4099077,Gluc_eco_rep1
        GSM4099078,Gluc_eco_rep2
        ...

    Flujo
    -----
    1. Lee el CSV de mapeo.
    2. Verifica que contenga las columnas indicadas.
    3. Reordena la matriz de conteos según el orden de `columna_matriz`.
    4. Renombra las columnas usando `columna_matriz` como llave
       y `columna_nuevo` como nombres nuevos.

    Parámetros
    ----------
    matriz_conteos : pandas.DataFrame
        Matriz de conteos (genes x muestras) con las columnas originales.
    ruta_mapeo : str o pathlib.Path
        Ruta al CSV con el mapeo de nombres de muestras.
    columna_matriz : str
        Nombre de la columna en el CSV que contiene los nombres actuales
        de las muestras (deben coincidir con las columnas de `matriz_conteos`).
    columna_nuevo : str
        Nombre de la columna en el CSV con los nombres nuevos que se
        asignarán a las muestras.

    Returns
    -------
    pandas.DataFrame
        Matriz de conteos con las columnas reordenadas y renombradas.

    Raises
    ------
    FileNotFoundError
        Si el archivo de mapeo no existe.
    KeyError
        Si las columnas de mapeo indicadas no están presentes en el CSV.
    """
    ruta_mapeo = Path(ruta_mapeo)
    if not ruta_mapeo.exists():
        raise FileNotFoundError(f"No se encontró el archivo de mapeo: {ruta_mapeo}")

    tabla_mapeo = pd.read_csv(ruta_mapeo)

    if columna_matriz not in tabla_mapeo.columns or columna_nuevo not in tabla_mapeo.columns:
        raise KeyError(
            f"Las columnas '{columna_matriz}' y/o '{columna_nuevo}' no están en el "
            f"archivo de mapeo. Columnas disponibles: {list(tabla_mapeo.columns)}"
        )

    nombres_en_mapeo = set(tabla_mapeo[columna_matriz])
    nombres_en_matriz = set(matriz_conteos.columns)

    nombres_faltantes = nombres_en_mapeo - nombres_en_matriz
    if nombres_faltantes:
        logger.error(
            "Algunas muestras del archivo de mapeo no se encuentran en la matriz "
            f"de conteos: {sorted(nombres_faltantes)}"
        )

    # Reordenar columnas según el orden de columna_matriz
    orden_original = list(tabla_mapeo[columna_matriz])
    columnas_presentes = [c for c in orden_original if c in matriz_conteos.columns]
    matriz_subconjunto = matriz_conteos.loc[:, columnas_presentes].copy()

    # Renombrar columnas
    diccionario_renombre = dict(zip(tabla_mapeo[columna_matriz], tabla_mapeo[columna_nuevo]))
    matriz_subconjunto = matriz_subconjunto.rename(columns=diccionario_renombre)

    logger.info(
        "Columnas de la matriz de conteos renombradas según archivo de mapeo "
        f"{ruta_mapeo}."
    )

    return matriz_subconjunto



# CPM y filtro por CPM

def cpm(matriz_conteos: pd.DataFrame):
    """
    Calcula los counts per million (CPM) por muestra para cada gen.

    Parámetros
    ----------
    matriz_conteos : pandas.DataFrame
        Matriz de conteos con genes en filas y muestras en columnas.

    Returns
    -------
    pandas.DataFrame
        Matriz de CPM con misma forma que `matriz_conteos`.
    """
    profundidad_por_muestra = matriz_conteos.sum(axis=0)
    profundidad_por_muestra = profundidad_por_muestra.replace(0, np.nan)

    matriz_cpm = (matriz_conteos.div(profundidad_por_muestra, axis=1)) * 1_000_000
    return matriz_cpm


def filtrar_por_cpm(matriz_conteos: pd.DataFrame, umbral_cpm: float, min_muestras: int, ruta_salida: str | Path | None = None,): 
    """
    Filtra genes según un umbral de CPM en al menos cierto número de muestras.

    Parámetros
    ----------
    matriz_conteos : pandas.DataFrame
        Matriz de conteos con genes en filas y muestras en columnas.
    umbral_cpm : float
        Umbral mínimo de CPM para considerar que un gen está "expresado"
        en una muestra.
    min_muestras : int
        Número mínimo de muestras en las que el CPM del gen debe ser
        mayor o igual a `umbral_cpm` para conservarlo.
    ruta_salida : str o pathlib.Path, opcional
        Si se proporciona, la matriz filtrada se guardará en esta ruta
        como archivo CSV.

    Returns
    -------
    pandas.DataFrame
        Matriz de conteos filtrada que contiene sólo los genes que
        cumplen el criterio de CPM.

    Notas
    -----
    Esta función es la que se usa en `main_deseq_fit.py` para aplicar el
    filtro previo al análisis de expresión diferencial.
    """
    matriz_cpm = cpm(matriz_conteos)

    num_muestras_con_cpm_suficiente = (matriz_cpm >= umbral_cpm).sum(axis=1)
    genes_pasan_filtro_cpm = num_muestras_con_cpm_suficiente >= min_muestras

    matriz_filtrada = matriz_conteos.loc[genes_pasan_filtro_cpm].copy()

    logger.info(
        f"Aplicando filtro CPM >= {umbral_cpm} en al menos {min_muestras} muestras: "
        f"{matriz_filtrada.shape[0]} genes restantes."
    )

    if ruta_salida is not None:
        ruta_salida = Path(ruta_salida)
        ruta_salida.parent.mkdir(parents=True, exist_ok=True)
        matriz_filtrada.to_csv(ruta_salida)
        logger.info(
            f"Matriz de conteos filtrada guardada en: {ruta_salida.resolve()}"
        )

    return matriz_filtrada


# Crear DeseqDataSet a partir de counts_T + metadata explícita

def crear_objeto_dds( matriz_conteos_transpuesta: pd.DataFrame, tabla_metadata: pd.DataFrame, nombre_columna_condicion: str = "condition",):
    """
    Crea un DeseqDataSet a partir de la matriz de conteos transpuesta y una
    tabla de metadata explícita (por ejemplo leída de metadata_ecoli.csv).

    Parámetros
    ----------
    matriz_conteos_transpuesta : pandas.DataFrame
        Matriz de conteos con muestras en filas y genes en columnas.
    tabla_metadata : pandas.DataFrame
        Tabla de metadata con muestras en el índice (mismo orden e índice
        que `matriz_conteos_transpuesta`).
    nombre_columna_condicion : str, opcional
        Nombre de la columna en la metadata que indica la condición
        experimental (por defecto 'condition').

    Returns
    -------
    pydeseq2.dds.DeseqDataSet
        Objeto listo para ejecutar el flujo DESeq2 (`dds.deseq2()`).

    Raises
    ------
    ValueError
        Si la matriz o la metadata están vacías, si sus índices no coinciden
        o si sólo se detecta una condición.
    """
    if matriz_conteos_transpuesta.empty:
        raise ValueError("La matriz de conteos transpuesta está vacía.")

    if tabla_metadata.empty:
        raise ValueError("La tabla de metadata está vacía.")

    if not matriz_conteos_transpuesta.index.equals(tabla_metadata.index):
        raise ValueError(
            "El índice de la matriz de conteos y la metadata no coincide. "
            "Ambas deben tener las mismas muestras en el mismo orden."
        )

    if nombre_columna_condicion not in tabla_metadata.columns:
        raise ValueError(
            f"La metadata no tiene la columna '{nombre_columna_condicion}'. "
            f"Columnas disponibles: {list(tabla_metadata.columns)}"
        )

    condiciones_unicas = tabla_metadata[nombre_columna_condicion].unique()
    if len(condiciones_unicas) < 2:
        raise ValueError(
            f"Solo se detectó una condición ({condiciones_unicas}). "
            "Se requieren al menos dos condiciones."
        )

    logger.info(
        f"Creando DeseqDataSet con {matriz_conteos_transpuesta.shape[0]} muestras, "
        f"{matriz_conteos_transpuesta.shape[1]} genes y condiciones: {condiciones_unicas}."
    )

    dds = DeseqDataSet(
        counts=matriz_conteos_transpuesta,
        metadata=tabla_metadata,
        design_factors=nombre_columna_condicion,
    )

    logger.info(
        f"DeseqDataSet creado: n_obs={dds.n_obs}, n_vars={dds.n_vars}."
    )

    return dds



#  Ejecutar DESeq2

def deseq2(dds: DeseqDataSet):
    """
    Ejecuta el análisis DESeq2 sobre un DeseqDataSet ya preparado.

    Parámetros
    ----------
    dds : pydeseq2.dds.DeseqDataSet
        Objeto creado con `crear_objeto_dds`.

    Returns
    -------
    pydeseq2.dds.DeseqDataSet
        El mismo objeto de entrada, pero con los resultados del modelo
        DESeq2 ya calculados (size factors, dispersions, etc.).

    Raises
    ------
    ValueError
        Si `dds` es None.
    """
    if dds is None:
        raise ValueError(
            "El objeto dds es None. Crea primero el objeto con crear_objeto_dds()."
        )

    logger.info("Iniciando análisis DESeq2 (dds.deseq2()).")
    dds.deseq2()
    logger.info(
        "Análisis DESeq2 completado. "
        f"n_obs (muestras) = {dds.n_obs}, n_vars (genes) = {dds.n_vars}."
    )

    if "size_factors" in dds.obs.columns:
        logger.info("Factores de tamaño (primeras filas):")
        cols = ["size_factors"]
        if "condition" in dds.obs.columns:
            cols.insert(0, "condition")
        logger.info(dds.obs[cols].head().to_string())

    if "dispersions" in dds.var.columns:
        logger.info("Dispersions (primeros genes):")
        logger.info(dds.var["dispersions"].head().to_string())

    return dds



# Estadísticas de DE (DeseqStats)

def estadisticas(dds: DeseqDataSet,condicion_prueba: str,condicion_control: str,nombre_columna_condicion: str = "condition",):
    """
    Calcula estadísticas de expresión diferencial para un contraste específico.

    Parámetros
    ----------
    dds : pydeseq2.dds.DeseqDataSet
        Objeto con el modelo DESeq2 ya ajustado (`deseq2(dds)`).
    condicion_prueba : str
        Nombre de la condición de prueba (tratamiento) tal como aparece
        en la columna `nombre_columna_condicion` de la metadata.
    condicion_control : str
        Nombre de la condición de referencia (control).
    nombre_columna_condicion : str, opcional
        Nombre de la columna de la metadata usada como factor de diseño
        (por defecto 'condition').

    Returns
    -------
    pandas.DataFrame
        DataFrame con las estadísticas del contraste (baseMean, log2FoldChange,lfcSE, stat, pvalue, padj, etc.).
    """
    logger.info(
        f"Iniciando análisis estadístico DESeq2: "
        f"{condicion_prueba} vs {condicion_control}."
    )

    ds = DeseqStats(
        dds,
        contrast=[nombre_columna_condicion, condicion_prueba, condicion_control],
    )
    ds.summary()

    tabla_resultados = ds.results_df.copy()

    logger.info("Análisis estadístico completado.")
    logger.info(f"Genes en resultados: {tabla_resultados.shape[0]}.")

    return tabla_resultados



# Clasificar genes en UP / DOWN / Non-DE

def clasificar_genes( res_df: pd.DataFrame,lfc_col: str,padj_col: str,padj_umbral: float,lfc_umbral: float,nueva_columna: str = "DE_status",):
    """
    Clasifica genes según log2FC y padj en tres categorías: 'UP', 'DOWN'
    y 'Non-DE'.

    Parámetros
    ----------
    res_df : pandas.DataFrame
        Tabla de resultados crudos de DESeq2 (salida de `estadisticas` /
        `main_deseq_fit.py`), que contiene al menos las columnas del
        log2FoldChange y del padj.
    lfc_col : str
        Nombre de la columna con el log2FoldChange.
    padj_col : str
        Nombre de la columna con el p-valor ajustado (padj).
    padj_umbral : float
        Umbral máximo de padj para considerar un gen como DE.
    lfc_umbral : float
        Umbral mínimo de la magnitud del log2FoldChange (|log2FC|) para
        considerar un gen como DE.
    nueva_columna : str, opcional
        Nombre de la columna nueva que contendrá la clasificación
        ('UP', 'DOWN' o 'Non-DE'). Por defecto 'DE_status'.

    Returns
    -------
    tabla_clasificada : pandas.DataFrame
        Copia de `res_df` con una columna adicional (`nueva_columna`) que
        indica la categoría de cada gen.
    resumen_categorias : dict
        Diccionario con el conteo de genes en cada categoría, por ejemplo:
        {'Non-DE': 1234, 'UP': 50, 'DOWN': 45}.

    Raises
    ------
    KeyError
        Si las columnas `lfc_col` o `padj_col` no se encuentran en `res_df`.
    """
    if lfc_col not in res_df.columns:
        raise KeyError(
            f"La columna '{lfc_col}' no se encuentra en los resultados."
        )
    if padj_col not in res_df.columns:
        raise KeyError(
            f"La columna '{padj_col}' no se encuentra en los resultados."
        )

    tabla = res_df.copy()

    tabla[lfc_col] = pd.to_numeric(
        tabla[lfc_col], errors="coerce"
    )
    tabla[padj_col] = pd.to_numeric(
        tabla[padj_col], errors="coerce"
    )

    condiciones_clasificacion = [
        (tabla[padj_col] <= padj_umbral) & (tabla[lfc_col] >= lfc_umbral),    # UP
        (tabla[padj_col] <= padj_umbral) & (tabla[lfc_col] <= -lfc_umbral),   # DOWN
    ]
    etiquetas = ["UP", "DOWN"]

    tabla[nueva_columna] = np.select(
        condiciones_clasificacion, etiquetas, default="Non-DE"
    )

    resumen_categorias = tabla[nueva_columna].value_counts().to_dict()
    logger.info(f"Resumen clasificación genes: {resumen_categorias}")

    return tabla, resumen_categorias



# 7. Guardar resultados

def guardar_resultados(df: pd.DataFrame,path_salida: str | Path,sep: str = ",",index: bool = True,):
    """
    Guarda un DataFrame de resultados en disco.

    Esta función se usa tanto en `main_deseq_fit.py` como en
    `main_deseq_clasificar.py`.

    Parámetros
    ----------
    df : pandas.DataFrame
        Tabla de resultados que se quiere guardar (por ejemplo, la tabla
        de DESeq2 o la tabla clasificada).
    path_salida : str o pathlib.Path
        Ruta del archivo de salida (se recomienda extensión .csv o .tsv).
    sep : str, opcional
        Separador de columnas a usar en el archivo (por defecto ',').
    index : bool, opcional
        Si es True, se guarda el índice del DataFrame en el archivo.
        Por defecto True.

    Returns
    -------
    None
        No retorna nada; sólo escribe el archivo en disco.
    """
    path_salida = Path(path_salida)
    path_salida.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path_salida, sep=sep, index=index)
    logger.info(f"Resultados guardados en: {path_salida.resolve()}")
