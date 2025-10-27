import argparse
from Bio import Entrez
import pandas as pd
from ftplib import FTP
from urllib.parse import urlparse
from pathlib import Path
from modules.descarga import assemblies_por_organismo, descargar_ensamblado

Entrez.email = "miryamzj@lcg.unam.mx"

def main():
    parser = argparse.ArgumentParser(
        description="Descarga archivos genómicos (fna, gff, faa) desde NCBI por organismo y accession ID."
    )
    parser.add_argument("--organismo", required=True, help="Nombre científico del organismo (entre comillas si tiene espacios).")
    parser.add_argument("--accession", required=True, help="ID de ensamblado (por ejemplo, GCF_000002945.2).")
    parser.add_argument("--outdir", default="descargas", help="Carpeta de salida donde guardar los archivos.")
    parser.add_argument("--incluir", nargs="+", default=["fna", "gff"], choices=["fna", "gff", "faa"],
                        help="Tipos de archivos a descargar (por defecto: fna gff).")
    parser.add_argument("--retmax", type=int, default=50, help="Número máximo de ensamblados a buscar (por defecto: 50).")

    args = parser.parse_args()

    print(f"\n Buscando ensamblados de {args.organismo}...")
    df = assemblies_por_organismo(args.organismo, retmax=args.retmax)

    print("\nColumnas disponibles:", df.columns.tolist())
    print(df.head(3).to_string())
    

    print(f"Descargando ensamblado {args.accession}...")
    descargar_ensamblado(df, args.accession, args.outdir, incluir=args.incluir)


if __name__ == "__main__":
    main()