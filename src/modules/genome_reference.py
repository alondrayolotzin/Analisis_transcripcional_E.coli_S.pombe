"""
genome_reference.py
Autora: Miryam Zamora
Fecha: 2025-10-28

Funciones auxiliares para buscar y descargar genomas de referencia desde NCBI.

Este módulo permite:

- Configurar el email para las consultas a NCBI (Entrez).
- Buscar ensamblados recientes para un organismo en la base de datos "assembly".
- Localizar una fila específica de ensamblado por accession.
- Listar archivos disponibles en el FTP de NCBI para un ensamblado.
- Filtrar archivos de interés (fna, gff, faa).
- Descargar los archivos seleccionados desde el FTP al disco local.

Está diseñado para ser utilizado por el script principal `main_download.py`
dentro de un pipeline de RNA-seq / genómica.
"""

from Bio import Entrez
import pandas as pd
from ftplib import FTP
from urllib.parse import urlparse
from pathlib import Path
import logging

# ========================
#  FUNCIONES BASE
# ========================
def set_email(email):
    Entrez.email = email
    
def assemblies_por_organismo(org, retmax=50):
    """
    Recupera un DataFrame con los ensamblados de referencia de un organismo.

    La búsqueda se hace en la base de datos "assembly", filtrando por:
    - Organismo (`org`).
    - Ensamblados marcados como "latest".
    - Ensamblados de referencia o representativos
      ("reference genome" OR "representative genome").

    Parameters
    ----------
    org : str
        Nombre del organismo (por ejemplo, "Escherichia coli").
    retmax : int, optional
        Número máximo de ensamblados a recuperar (por defecto: 50).

    Returns
    -------
    pandas.DataFrame
        DataFrame con una fila por ensamblado y las columnas retornadas por
        `Entrez.esummary` (DocumentSummary). Si no se encuentran ensamblados,
        el DataFrame estará vacío.
    """
    
    term = f'({org}[Organism]) AND latest[filter] AND ("reference genome"[filter] OR "representative genome"[filter])'
    with Entrez.esearch(db="assembly", term=term, retmax=retmax) as h:
        res = Entrez.read(h)
    uids = res.get("IdList", [])

    filas = []
    for uid in uids:
        with Entrez.esummary(db="assembly", id=uid) as h:
            s = Entrez.read(h)
        filas.append(s["DocumentSummarySet"]["DocumentSummary"][0])

    df = pd.DataFrame(filas)
    logging.info(f"Se encontraron {len(df)} ensamblados para {org}")
    return df


def _row_por_accession(df, accession):
    """
    Devuelve la fila del ensamblado que coincide con el accession solicitado.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame con información de ensamblados (salida de assemblies_por_organismo).
    accession : str
        Accession del ensamblado (por ejemplo, "GCF_000005845.2").

    Returns
    -------
    pandas.Series
        Fila del DataFrame correspondiente al ensamblado solicitado.

    Raises
    ------
    ValueError
        Si el accession no se encuentra en el DataFrame.
    """
    
    mask = df["AssemblyAccession"].eq(accession)
    if not mask.any():
        logging.error(f"No se encontró {accession} en el DataFrame de ensamblados.")
        raise ValueError(f"No se encontró {accession} en el DataFrame.")
    return df.loc[mask].iloc[0]


def _ftp_listar(ftp_url):
    """
    Lista los archivos disponibles en el FTP de NCBI para un ensamblado.

    Parameters
    ----------
    ftp_url : str
        URL FTP al directorio del ensamblado en NCBI, por ejemplo:
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/...".

    Returns
    -------
    list of str
        Lista de nombres de archivo presentes en ese directorio FTP.
    """
    p = urlparse(ftp_url)
    with FTP(p.hostname, timeout=60) as ftp:
        ftp.login()
        ftp.cwd(p.path)
        return ftp.nlst()


def _targets(files, incluir=("fna", "gff", "faa")):
    """
    Filtra los archivos de un ensamblado según los sufijos de interés.

    Parameters
    ----------
    files : list of str
        Lista de nombres de archivo en el directorio FTP del ensamblado.
    incluir : tuple of str, optional
        Tipos de archivo a incluir. Valores posibles:
        - "fna" → archivos con sufijo "_genomic.fna.gz"
        - "gff" → archivos con sufijo "_genomic.gff.gz"
        - "faa" → archivos con sufijo "_protein.faa.gz"

    Returns
    -------
    list of str
        Subconjunto de `files` que termina con los sufijos seleccionados.
    """
    suf = {
        "fna": "_genomic.fna.gz",
        "gff": "_genomic.gff.gz",
        "faa": "_protein.faa.gz"
    }
    return [f for f in files for t in incluir if f.endswith(suf[t])]


def _ftp_descargar(ftp_url, nombres, outdir):
    """
    Descarga los archivos seleccionados desde el FTP de NCBI a un directorio local.

    Parameters
    ----------
    ftp_url : str
        URL FTP al directorio del ensamblado en NCBI.
    nombres : list of str
        Lista de nombres de archivo a descargar desde ese directorio.
    outdir : str or pathlib.Path
        Directorio local donde se guardarán los archivos descargados.
        Se creará si no existe.
    """

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

    logging.info(f"Descargas completadas en: {outdir.resolve()}")

def descargar_ensamblado(df, accession, outdir, incluir=("fna", "gff", "faa")):
    """
    Descarga los archivos genómicos asociados a un ensamblado concreto.

    A partir de un DataFrame con información de ensamblados (por ejemplo, la
    salida de `assemblies_por_organismo`), localiza la fila correspondiente
    al accession solicitado, identifica la ruta FTP en RefSeq o GenBank y
    descarga los archivos seleccionados (fna/gff/faa).

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame con información de ensamblados ("assembly").
    accession : str
        Accession del ensamblado a descargar (por ejemplo "GCF_000005845.2").
    outdir : str or pathlib.Path
        Directorio local donde guardar los archivos descargados.
    incluir : tuple of str, optional
        Tipos de archivos a descargar. Valores posibles: "fna", "gff", "faa".
        Por defecto, ("fna", "gff", "faa").

    Returns
    -------
    (pandas.Series, str, list of str)
        Tupla con:
        - Fila del DataFrame correspondiente al ensamblado.
        - URL FTP utilizada.
        - Lista de nombres de archivo descargados.

    Raises
    ------
    ValueError
        Si el accession no se encuentra en el DataFrame.
    """
    fila = _row_por_accession(df, accession)
    ftp_url = fila.get("FtpPath_RefSeq") or fila.get("FtpPath_GenBank")

    org_name = fila.get('Organism') or fila.get('OrganismName')
    acc_name = fila.get('AssemblyAccession') or fila.get('AssemblyAccn') or fila.get('Accession')
    logging.info(f"Iniciando descarga de {org_name} ({acc_name})")

    files = _ftp_listar(ftp_url)
    names = _targets(files, incluir=incluir)
    _ftp_descargar(ftp_url, names, outdir)

    logging.info(f"Descarga completada: {org_name} ({acc_name})")
    return fila, ftp_url, names 
