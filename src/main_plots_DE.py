"""
main_plots_DE.py
Autora: Alondra Márquez
Fecha: 2025-11-24

Descripción general
-------------------
Script de línea de comandos para generar las gráficas posteriores al análisis
de expresión diferencial a partir de archivos ya procesados.

Este script:

Lee la tabla de resultados clasificados de DESeq2
    (por ejemplo: deseq2_resultados_clasificados.csv).
Genera un volcano plot para el contraste especificado
    (condición de prueba vs condición control).
Genera un barplot con el número de genes UP / DOWN / Non-DE.
Genera un heatmap de los genes más diferencialmente expresados, integrando:
    - Matriz de conteos filtrados por CPM.
    - Metadata de las muestras (sample_id, condition, etc.).

Entradas esperadas
------------------
Archivo de resultados DESeq2 clasificados:
Debe contener al menos las columnas: 'gene_id', 'log2FoldChange', 'padj', 'DE_status'.
Matriz de conteos filtrados por CPM:
    - Formato genes x muestras (CSV).
    - Puede incluir una columna 'gene_id' que será usada como índice.
Metadata de las muestras:
    - CSV con columnas 'sample_id' y 'condition' (mínimo).
    - 'sample_id' debe coincidir con los nombres de las columnas de la matriz de conteos.

Uso 
----------
    python main_plots_DE.py \
        --organismo "E. coli" \
        --resultados results/E_coli/deseq2_resultados_clasificados.csv \
        --counts results/E_coli/counts_filtrados_cpm.csv \
        --metadata data/E_coli/metadata_ecoli.csv \
        --condicion-prueba "Glicerol" \
        --condicion-control "Glucosa" \
        --outdir results/E_coli/figuras \
        --padj-umbral 0.05 \
        --lfc-umbral 1.0 \
        --top-genes-heatmap 50

"""
import argparse
from pathlib import Path

from modules.plots_DE import (
    leer_resultados_deseq,
    hacer_volcano,
    hacer_barplot_up_down,
    hacer_heatmap_top_genes,
)


def main():
    parser = argparse.ArgumentParser(description=( "Generación de gráficas post-DESeq2 (volcano, barplot, heatmap") )

    parser.add_argument("--organismo", required=True,help="Nombre del organismo (ej. 'E. coli', 'S. pombe').",)
    parser.add_argument("--resultados",required=True,help="Ruta a deseq2_resultados_clasificados.csv.",)
    parser.add_argument("--counts",required=True,help="Ruta a counts_filtrados_cpm.csv (genes x muestras).",)
    parser.add_argument("--metadata",required=True,help="Ruta a metadata CSV (con columnas 'sample_id' y 'condition').",)
    parser.add_argument("--condicion-prueba",required=True,help="Nombre de la condición de prueba (ej. 'Glicerol').",)
    parser.add_argument("--condicion-control",required=True,help="Nombre de la condición control (ej. 'Glucosa').",)
    parser.add_argument("--outdir",required=True,help="Directorio de salida para las figuras.",)
    parser.add_argument("--padj-umbral",type=float,default=0.05,help="Umbral de padj para líneas de referencia en volcano (por defecto 0.05).",)
    parser.add_argument("--lfc-umbral",type=float, default=1.0, help="Umbral de |log2FC| para líneas de referencia en volcano (por defecto 1.0).",)
    parser.add_argument("--top-genes-heatmap",type=int,default=100,help="Número de genes top DE para el heatmap (por defecto 50).",)

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Leer resultados clasificados
    res_df = leer_resultados_deseq(args.resultados, sep=",")

    # 2) Volcano plot
    volcano_path = outdir / f"volcano_{args.organismo.replace(' ', '_')}.png"
    hacer_volcano(
        res_df=res_df,
        organismo=args.organismo,
        condicion_prueba=args.condicion_prueba,
        condicion_control=args.condicion_control,
        padj_umbral=args.padj_umbral,
        lfc_umbral=args.lfc_umbral,
        outpath=volcano_path,
    )

    # 3) Barplot UP/DOWN/Non-DE
    barplot_path = outdir / f"barplot_DE_{args.organismo.replace(' ', '_')}.png"
    hacer_barplot_up_down(
        res_df=res_df,
        organismo=args.organismo,
        outpath=barplot_path,
    )
    # 4) Heatmap de top genes DE
    heatmap_path = outdir / f"heatmap_top{args.top_genes_heatmap}_{args.organismo.replace(' ', '_')}.png"
    hacer_heatmap_top_genes(
        counts_path=args.counts,
        resultados_df=res_df,
        metadata_path=args.metadata,
        organismo=args.organismo,
        condicion_col="condition",
        status_col="DE_status",
        n_genes=args.top_genes_heatmap,
        sep_counts=",",
        sep_metadata=",",
        outpath=heatmap_path,
    )


if __name__ == "__main__":
    main()
