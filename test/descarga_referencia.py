from Bio import Entrez
import pandas as pd
from io import StringIO

Entrez.email = "miryamzj@lcg.unam.mx"

def assemblies_por_organismo(org, retmax=50):
    term = f'({org}[Organism]) AND latest[filter] AND ("reference genome"[filter] OR "representative genome"[filter])'
    with Entrez.esearch(db="assembly", term=term, retmax=retmax) as h:
        res = Entrez.read(h)
    uids = res["IdList"]

    filas = []
    for uid in uids:
        with Entrez.esummary(db="assembly", id=uid) as h:
            s = Entrez.read(h)
        filas.append(s["DocumentSummarySet"]["DocumentSummary"][0])

    return pd.DataFrame(filas)

# 1) DEFINE los nombres de cada organismo (sin el strain en [Organism])
S_POMBE = "Schizosaccharomyces pombe"
E_COLI  = "Escherichia coli"

# 2) LLAMA la función con cada nombre
df_pombe = assemblies_por_organismo(S_POMBE, retmax=20)
df_coli  = assemblies_por_organismo(E_COLI,  retmax=50)

# 3) Muestra resultados
print("Pombe:", len(df_pombe), "ensamblajes")
print(df_pombe[["Organism","AssemblyName","AssemblyAccession","RefSeq_category"]].head(10).to_string(index=False))

print("\nE. coli:", len(df_coli), "ensamblajes")
print(df_coli[["Organism","AssemblyName","AssemblyAccession","RefSeq_category"]].head(10).to_string(index=False))


# ----------------------------------------

from ftplib import FTP
from urllib.parse import urlparse
from pathlib import Path

# --- helpers que ya usabas ---
def _row_por_accession(df, accession):
    mask = df["AssemblyAccession"].eq(accession)
    if not mask.any():
        raise ValueError(f"No se encontró {accession} en el DataFrame.")
    return df.loc[mask].iloc[0]

def _ftp_listar(ftp_url):
    p = urlparse(ftp_url)
    with FTP(p.hostname, timeout=60) as ftp:
        ftp.login()
        ftp.cwd(p.path)
        return ftp.nlst()

def _targets(files, incluir=("fna", "gff", "faa")):
    suf = {"fna":"_genomic.fna.gz", "gff":"_genomic.gff.gz", "faa":"_protein.faa.gz"}
    return [f for f in files for t in incluir if f.endswith(suf[t])]

def _ftp_descargar(ftp_url, nombres, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    p = urlparse(ftp_url)
    with FTP(p.hostname, timeout=60) as ftp:
        ftp.login()
        ftp.cwd(p.path)
        for name in nombres:
            dest = outdir / name
            print("↓", name, "→", dest)
            with open(dest, "wb") as fh:
                ftp.retrbinary(f"RETR {name}", fh.write)
    print("✔ Descargas en:", outdir.resolve())

# --- funciones específicas que ya tenías (sin prefer cases) ---
def descargar_ecoli_ref(df_coli, outdir, incluir=("fna","gff")):
    """
    E. coli K-12 MG1655 (reference): GCF_000005845.2
    """
    fila = _row_por_accession(df_coli, "GCF_000005845.2")
    ftp_url = fila.get("FtpPath_RefSeq") or fila.get("FtpPath_GenBank")
    print(f"\nE. coli → {fila['AssemblyName']} ({fila['AssemblyAccession']})")
    files  = _ftp_listar(ftp_url)
    names  = _targets(files, incluir=incluir)
    _ftp_descargar(ftp_url, names, outdir)
    return fila, ftp_url, names

def descargar_spombe_ref(df_pombe, outdir, incluir=("fna","gff")):
    """
    S. pombe referencia: GCF_000002945.2 (ASM294v3)
    """
    fila = _row_por_accession(df_pombe, "GCF_000002945.2")
    ftp_url = fila.get("FtpPath_RefSeq") or fila.get("FtpPath_GenBank")
    print(f"\nS. pombe → {fila['AssemblyName']} ({fila['AssemblyAccession']})")
    files  = _ftp_listar(ftp_url)
    names  = _targets(files, incluir=incluir)
    _ftp_descargar(ftp_url, names, outdir)
    return fila, ftp_url, names



# ESTO ES YA DESDE MAIN 

# /export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/src/main_download.py
from pathlib import Path
from Bio import Entrez
from modules.bioproject import (
    descargar_ecoli_ref, descargar_spombe_ref
)
from modules.tus_busquedas import assemblies_por_organismo  # <-- importa tu función donde la tengas

# Configura correo NCBI
Entrez.email = "miryamzj@lcg.unam.mx"

# Rutas: src → proyecto → data
SRC_DIR   = Path(__file__).resolve().parent
ROOT_DIR  = SRC_DIR.parent
DATA_DIR  = ROOT_DIR / "data" #Calcula la ruta al directorio data, ya no le tienes que copiar la ruta.
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Organismos
S_POMBE = "Schizosaccharomyces pombe"
E_COLI  = "Escherichia coli"

# Busca ensamblajes
df_pombe = assemblies_por_organismo(S_POMBE, retmax=20)
df_coli  = assemblies_por_organismo(E_COLI,  retmax=50)

# Subcarpetas de salida dentro de data/
OUT_POMBE = DATA_DIR / "S_pombe_ASM294v3"  #Aquí se supone que solito calcula la ruta a data 
OUT_COLI  = DATA_DIR / "E_coli_MG1655"     

# Descargas (FASTA + GFF por defecto)
descargar_spombe_ref(df_pombe, outdir=OUT_POMBE)
descargar_ecoli_ref(df_coli,  outdir=OUT_COLI)

print("\nTodo listo ✅")
print("Pombe →", OUT_POMBE.resolve())
print("E. coli →", OUT_COLI.resolve())
