from Bio import Entrez
import pandas as pd
from ftplib import FTP
from urllib.parse import urlparse
from pathlib import Path

Entrez.email = "miryamzj@lcg.unam.mx"

# ========================
#  FUNCIONES BASE
# ========================

def assemblies_por_organismo(org, retmax=50):
    """Devuelve un DataFrame con los ensamblados más recientes del organismo indicado."""
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


def _row_por_accession(df, accession):
    """Devuelve la fila del ensamblado que coincide con el accession."""
    mask = df["AssemblyAccession"].eq(accession)
    if not mask.any():
        raise ValueError(f"No se encontró {accession} en el DataFrame.")
    return df.loc[mask].iloc[0]


def _ftp_listar(ftp_url):
    """Lista los archivos disponibles en el FTP de NCBI para un ensamblado."""
    p = urlparse(ftp_url)
    with FTP(p.hostname, timeout=60) as ftp:
        ftp.login()
        ftp.cwd(p.path)
        return ftp.nlst()


def _targets(files, incluir=("fna", "gff", "faa")):
    """Filtra los archivos según los sufijos de interés."""
    suf = {
        "fna": "_genomic.fna.gz",
        "gff": "_genomic.gff.gz",
        "faa": "_protein.faa.gz"
    }
    return [f for f in files for t in incluir if f.endswith(suf[t])]


def _ftp_descargar(ftp_url, nombres, outdir):
    """Descarga los archivos seleccionados desde el FTP."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    p = urlparse(ftp_url)

    with FTP(p.hostname, timeout=60) as ftp:
        ftp.login()
        ftp.cwd(p.path)
        for name in nombres:
            dest = outdir / name
            print(f"Descargando {name} → {dest}")
            with open(dest, "wb") as fh:
                ftp.retrbinary(f"RETR {name}", fh.write)

    print(f"\n Descargas completadas en: {outdir.resolve()}\n")

def descargar_ensamblado(df, accession, outdir, incluir=("fna", "gff", "faa")):
    """
    Descarga los archivos genómicos asociados a un ensamblado (por accession ID).
    """
    fila = _row_por_accession(df, accession)
    ftp_url = fila.get("FtpPath_RefSeq") or fila.get("FtpPath_GenBank")

    print(f"\n Descargando {fila.get('Organism') or fila.get('OrganismName')} ({fila['AssemblyAccession'] or fila.get('AssemblyAccn') or fila.get('Accession')})")
    files = _ftp_listar(ftp_url)
    names = _targets(files, incluir=incluir)
    _ftp_descargar(ftp_url, names, outdir)

    return fila, ftp_url, names    