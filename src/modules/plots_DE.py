"""
plots_DE.py
Autor: Alondra Márquez
Fecha: 2025-11-24

Descripción general
-------------------
Este módulo agrupa las funciones encargadas de generar las gráficas
posteriores al análisis de expresión diferencial (DESeq2 / PyDESeq2).
Incluye rutinas para:

    • Leer la tabla de resultados clasificados de DESeq2
      (incluyendo log2FC, padj y DE_status).
    • Generar volcano plots para visualizar la magnitud del cambio de expresión
      frente a la significancia estadística.
    • Generar barplots con el número de genes UP / DOWN / Non-DE.
    • Generar heatmaps/clustermaps de los genes más diferencialmente expresados,
      integrando la matriz de conteos y la metadata de las muestras.

Este módulo es utilizado directamente por los scripts:

    - main_plots_DE.py (o scripts equivalentes)
        Encargados de leer las tablas de resultados y matrices de conteos,
        y llamar a estas funciones para producir las figuras finales.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ------------------------------------------------------------------
# Ajustes estéticos globales
# ------------------------------------------------------------------
# Tema general: fondo clarito, ejes limpios, tipografía un poco más grande
sns.set_theme(style="whitegrid", context="talk")

# Paletas de colores consistentes en todo el módulo
COLOR_UP = "#d73027"      # rojo ladrillo
COLOR_DOWN = "#4575b4"    # azul frío
COLOR_NONDE = "#bdbdbd"   # gris neutro

BAR_PALETTE = {
    "UP": COLOR_UP,
    "DOWN": COLOR_DOWN,
    "Non-DE": COLOR_NONDE,
}


def leer_resultados_deseq(path, sep=","):
    """
    Lee resultados clasificados de DESeq2 y verifica que tengan las columnas mínimas
    necesarias para generar las gráficas.
    """
    path = Path(path)
    df = pd.read_csv(path, sep=sep)

    # Eliminar índice viejo si existe
    if "Unnamed: 0" in df.columns:
        df = df.drop(columns=["Unnamed: 0"])

    # Si gene_id está como índice, rescatarlo
    if "gene_id" not in df.columns and df.index.name == "gene_id":
        df = df.reset_index()

    cols_necesarias = {"gene_id", "log2FoldChange", "padj", "DE_status"}
    faltan = cols_necesarias - set(df.columns)
    if faltan:
        raise ValueError(f"Faltan columnas en resultados DESeq2: {faltan}")

    return df


def hacer_volcano(
    res_df,
    organismo,
    condicion_prueba,
    condicion_control,
    lfc_col="log2FoldChange",
    padj_col="padj",
    status_col="DE_status",
    padj_umbral=0.05,
    lfc_umbral=1.0,
    outpath=None,
):
    """
    Genera un volcano plot a partir de resultados DESeq2 clasificados.
    """
    df = res_df.copy()

    # Evitar problemas con padj = 0 o NaN
    df = df.dropna(subset=[padj_col, lfc_col])
    df["minus_log10_padj"] = -np.log10(df[padj_col].replace(0, np.nan))

    plt.figure(figsize=(8, 6))

    # Colores por estado (paleta "bonita")
    palette = {
        "UP": COLOR_UP,
        "DOWN": COLOR_DOWN,
        "Non-DE": COLOR_NONDE,
    }

    # Dibujar primero Non-DE para que queden "al fondo"
    orden_status = ["Non-DE", "DOWN", "UP"]
    for status in orden_status:
        sub = df[df[status_col] == status]
        if sub.empty:
            continue

        plt.scatter(
            sub[lfc_col],
            sub["minus_log10_padj"],
            s=18 if status == "Non-DE" else 22,
            alpha=0.35 if status == "Non-DE" else 0.75,
            label=status,
            c=palette.get(status, COLOR_NONDE),
            edgecolors="none",
        )

    # Líneas de corte (un poco más discretas)
    corte_color = "#666666"
    plt.axvline(x=lfc_umbral, color=corte_color, linestyle="--", linewidth=1)
    plt.axvline(x=-lfc_umbral, color=corte_color, linestyle="--", linewidth=1)
    plt.axhline(y=-np.log10(padj_umbral), color=corte_color, linestyle="--", linewidth=1)

    plt.xlabel("log2(Fold Change)")
    plt.ylabel("-log10(padj)")
    plt.title(
        f"Volcano plot - {organismo}\n{condicion_prueba} vs {condicion_control}"
    )

    # Colocar la leyenda fuera para que no tape puntos
    plt.legend(
        title=status_col,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )

    # Limpiar bordes superior y derecho (más minimalista)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if outpath is not None:
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
        print(f"[OK] Volcano plot guardado en: {outpath}")

    plt.close()


def hacer_barplot_up_down(res_df, organismo, outpath=None, status_col="DE_status"):
    """
    Grafica el número de genes en cada categoría de expresión diferencial
    (UP, DOWN y Non-DE) usando un barplot simple.
    """
    conteo = res_df[status_col].value_counts().reindex(["UP", "DOWN", "Non-DE"])

    plt.figure(figsize=(5, 4))

    # Usar la paleta definida arriba
    colores = [BAR_PALETTE.get(cat, COLOR_NONDE) for cat in conteo.index]

    sns.barplot(
        x=conteo.index,
        y=conteo.values,
        palette=colores,
        edgecolor="black",
    )

    plt.ylabel("Número de genes")
    plt.xlabel("Categoría DE")
    plt.title(f"Resumen DE - {organismo}")

    # Mostrar valores encima de cada barra
    for i, v in enumerate(conteo.values):
        plt.text(
            i,
            v + max(conteo.values) * 0.01,
            str(v),
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if outpath is not None:
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
        print(f" Barplot UP/DOWN guardado en: {outpath}")

    plt.close()


def hacer_heatmap_top_genes(
    counts_path,
    resultados_df,
    metadata_path,
    organismo,
    condicion_col="condition",
    status_col="DE_status",
    n_genes=50,
    sep_counts=",",
    sep_metadata=",",
    outpath=None,
):
    """
    Genera un heatmap (clustermap) de los genes más diferencialmente expresados.
    """
    counts_path = Path(counts_path)
    metadata_path = Path(metadata_path)

    # 1) Leer matriz de conteos normalizados
    df_counts = pd.read_csv(counts_path, sep=sep_counts)
    if "gene_id" in df_counts.columns:
        df_counts = df_counts.set_index("gene_id")
    normCounts = df_counts

    # 2) Seleccionar genes significativos (UP y DOWN)
    df_de = resultados_df.copy()
    df_de = df_de.dropna(subset=["padj"])
    df_de = df_de[df_de[status_col].isin(["UP", "DOWN"])]

    df_de = df_de.sort_values("padj").head(n_genes)
    genesSig = df_de["gene_id"].unique()
    print(f"Genes seleccionados para heatmap: {len(genesSig)}")

    genes_en_matriz = normCounts.index.intersection(genesSig)
    if genes_en_matriz.empty:
        raise ValueError("Ninguno de los genes seleccionados está en la matriz de conteos.")

    heatmap_data = np.log1p(normCounts.loc[genes_en_matriz])
    print(f" Dimensiones de heatmap_data: {heatmap_data.shape}")

    # 4) Leer metadata y alinear muestras
    meta = pd.read_csv(metadata_path, sep=sep_metadata)
    if "sample_id" not in meta.columns:
        raise ValueError("La metadata debe tener una columna 'sample_id'.")

    meta = meta.set_index("sample_id")

    muestras_comunes = [m for m in heatmap_data.columns if m in meta.index]
    heatmap_data = heatmap_data[muestras_comunes]
    meta = meta.loc[muestras_comunes]

    # 5) Colores por condición
    condiciones_unicas = meta[condicion_col].unique()
    cond_palette = dict(
        zip(
            condiciones_unicas,
            sns.color_palette("Set2", n_colors=len(condiciones_unicas)),
        )
    )
    col_colors = meta[condicion_col].map(cond_palette)

    # 6) Clustermap con cmap diverging bonito
    g = sns.clustermap(
        heatmap_data,
        z_score=0,                # normalizar por gen
        cmap="RdBu_r",            # azul-bajo / rojo-alto
        col_cluster=True,
        row_cluster=True,
        col_colors=col_colors,
        figsize=(10, 8),
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),  # posición de la barra de color
    )

    # Título general
    g.fig.suptitle(f"Heatmap top {n_genes} genes DE - {organismo}", y=1.02)

    # Quitar el marco alrededor de la matriz para un look más limpio
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("gene_id")

    if outpath is not None:
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        g.savefig(outpath, dpi=300, bbox_inches="tight")
        print(f"Heatmap guardado en: {outpath}")

    plt.close()
