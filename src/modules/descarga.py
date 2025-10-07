from Bio import Entrez
import pandas as pd
from ftplib import FTP
from urllib.parse import urlparse
from pathlib import Path 

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

# --- helpers internos ---
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

# --- funciones específicas para descargar los genomas ---
def descargar_ecoli_ref(df_coli, outdir, incluir=("fna","gff")):
    fila = _row_por_accession(df_coli, "GCF_000005845.2")
    ftp_url = fila.get("FtpPath_RefSeq") or fila.get("FtpPath_GenBank")
    print(f"\nE. coli → {fila['AssemblyName']} ({fila['AssemblyAccession']})")
    files  = _ftp_listar(ftp_url)
    names  = _targets(files, incluir=incluir)
    _ftp_descargar(ftp_url, names, outdir)
    return fila, ftp_url, names

def descargar_spombe_ref(df_pombe, outdir, incluir=("fna","gff")):
    fila = _row_por_accession(df_pombe, "GCF_000002945.2")
    ftp_url = fila.get("FtpPath_RefSeq") or fila.get("FtpPath_GenBank")
    print(f"\nS. pombe → {fila['AssemblyName']} ({fila['AssemblyAccession']})")
    files  = _ftp_listar(ftp_url)
    names  = _targets(files, incluir=incluir)
    _ftp_descargar(ftp_url, names, outdir)
    return fila, ftp_url, names

