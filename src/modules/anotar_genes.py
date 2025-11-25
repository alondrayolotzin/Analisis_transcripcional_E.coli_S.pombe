"""
anotacion_funcional.py

Módulo para:

1. Leer una LISTA de genes ya filtrados (por ejemplo, genes_UP.txt o genes_DOWN.txt),
   uno por línea.
2. Convertir esos IDs a IDs de KEGG (eco:b0001, spo:XXXX, etc.).
   - Si los IDs ya son locus_tag compatibles con KEGG (b0001, b0002...), usa
     --source-db direct y solo se antepone el código de organismo.
   - Si vienen de otra BD (ncbi-geneid, uniprot, ...), se usa kegg_conv.
3. Consultar KEGG en CHUNKS usando kegg_get para obtener anotación básica:
   NAME, DEFINITION, KO, PATHWAYS.
4. Guardar un archivo TSV con:
   - gene_id        (ID limpio del archivo de entrada)
   - source_id      (ID usado para la conversión)
   - kegg_id        (eco:b0001, etc.)
   - name
   - definition
   - ko
   - pathways
   - set_label      (opcional: etiqueta UP / DOWN o la que des)

Ejemplo esperado de lista de genes (genes_UP.txt):

b0001
b0002
b0003
...

Ejemplo de uso:

nohup python src/modules/anotacion_funcional.py \
  --genes-file results/genes_UP.txt \
  --org-code eco \
  --source-db direct \
  --set-label UP \
  --out results/eco_genes_UP_annotation.tsv \
  > eco_genes_UP_annotation.out 2>&1 &
"""

import argparse
import logging
import time
from typing import Dict, List
import re
from pathlib import Path
import pandas as pd
from Bio.KEGG import REST


# ---------------------------------------------------------------------
# Utilidades para leer y limpiar IDs
# ---------------------------------------------------------------------


def cargar_lista_genes(path: str) -> List[str]:
    """
    Lee una lista de genes desde un archivo de texto, uno por línea.

    Ignora líneas vacías y líneas que empiezan con '#'.
    """
    genes: List[str] = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            genes.append(s)
    logging.info("Lista de genes cargada desde %s (n = %d)", path, len(genes))
    return genes


def limpiar_id_gen(gene_id: str) -> str:
    """
    Limpia el ID del gen para quitar prefijos tipo 'ID=', 'cds-' o 'SPOM_'.

    También normaliza la letra final 'C' -> 'c' en IDs de S. pombe,
    por ejemplo: SPOM_SPAC1002.07C -> SPAC1002.07c.
    Para E. coli, b0001 se queda como b0001.
    """
    gid = gene_id.strip()

    # Cosas generales
    if gid.startswith("ID="):
        gid = gid.split("=", 1)[1]
    if gid.startswith("cds-"):
        gid = gid[4:]

    # Prefijo específico de S. pombe
    if gid.startswith("SPOM_"):
        # SPOM_SPAC1002.07C -> SPAC1002.07C
        gid = gid.split("SPOM_", 1)[1]

    # Si termina en letra después de un número (ej. ...07C), bajar a minúscula
    if len(gid) >= 2 and gid[-1].isalpha() and gid[-2].isdigit():
        gid = gid[:-1] + gid[-1].lower()

    return gid



# ---------------------------------------------------------------------
# Conversión de IDs a IDs de KEGG
# ---------------------------------------------------------------------
_SPO_MAP_CACHE: Dict[str, str] | None = None

_SPO_MAP_CACHE: Dict[str, str] | None = None

def cargar_mapa_spombe(path: str = "data/spo_kegg_genes.list") -> Dict[str, str]:
    """
    Lee el archivo spo_kegg_genes.list y construye un mapa:

        ID_normalizado (SPACxxxx.xx[c], SPBCxx..., etc.) -> spo:NUM

    Se cachea en memoria para no leerlo muchas veces.
    """
    global _SPO_MAP_CACHE
    if _SPO_MAP_CACHE is not None:
        return _SPO_MAP_CACHE

    mapa: Dict[str, str] = {}

    ruta = Path(path)
    if not ruta.exists():
        raise FileNotFoundError(
            f"No se encontró el archivo de mapa de S. pombe: {ruta} "
            "(ajusta la ruta en cargar_mapa_spombe)."
        )

    with ruta.open() as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            kegg_id = cols[0].strip()      # ej. 'spo:2541932'
            resto = " ".join(cols[1:])     # el resto de la línea

            # Partimos el resto en “tokens” separados por espacios, comas, etc.
            # Buscamos cosas que parezcan IDs de PomBase (contienen 'SP' y '.')
            tokens = (
                resto.replace(",", " ")
                     .replace(";", " ")
                     .replace("/", " ")
                     .split()
            )

            for tok in tokens:
                if "SP" in tok and "." in tok:
                    # Ejemplos de tok:
                    #   SPOM_SPAC212.11
                    #   SPAC11D3.03C
                    #   SPBC19G7.04c
                    alias_norm = limpiar_id_gen(tok)
                    mapa[alias_norm] = kegg_id

    _SPO_MAP_CACHE = mapa
    logging.info("Mapa S. pombe cargado: %d IDs", len(mapa))
    return mapa


def convertir_ids_a_kegg(
    gene_ids: List[str],
    source_db: str,
    org_code: str,
    delay: float = 0.5,
) -> pd.DataFrame:
    """
    Convierte IDs de genes a IDs de KEGG.

    Modos:

    - source_db = 'direct' / 'locus_tag' / 'kegg':
      Los IDs de la lista ya son los que KEGG usa como parte derecha
      del ID (por ejemplo, b0001 -> eco:b0001). No se llama a kegg_conv.

    - source_db = algo como 'ncbi-geneid', 'ncbi-proteinid', 'uniprot', etc.:
      Se llama a REST.kegg_conv(org_code, source_id) por cada gen.
    """
    registros: List[Dict[str, str]] = []

    # -------- CASO 1: IDs directos / locus_tag --------
    if source_db.lower() in {"direct", "locus_tag", "kegg"}:
        logging.info(
            "Usando IDs directos (%s) para org_code=%s",
            source_db,
            org_code,
        )

        # Caso especial: S. pombe (spo) -> necesitamos mapa SPAC -> spo:NUM
        if org_code == "spo":
            mapa_spo = cargar_mapa_spombe()
            for raw_id in gene_ids:
                clean_id = limpiar_id_gen(raw_id)
                if not clean_id:
                    continue

                kegg_id = mapa_spo.get(clean_id)
                if not kegg_id:
                    logging.warning(
                        "No se encontró ID KEGG para S. pombe: %s",
                        clean_id,
                    )
                    continue

                registros.append(
                    {
                        "gene_id": clean_id,
                        "source_id": clean_id,
                        "kegg_id": kegg_id,
                    }
                )

        else:
            # Resto de organismos (p.ej. E. coli con b0001 -> eco:b0001)
            for raw_id in gene_ids:
                clean_id = limpiar_id_gen(raw_id)
                if not clean_id:
                    continue

                kegg_id = f"{org_code}:{clean_id}"
                registros.append(
                    {
                        "gene_id": clean_id,
                        "source_id": clean_id,
                        "kegg_id": kegg_id,
                    }
                )

        df_ids = pd.DataFrame(registros)
        logging.info(
            "Conversión 'directa' a IDs KEGG completada para %d genes",
            df_ids.shape[0],
        )
        return df_ids

    # -------- CASO 2: usar kegg_conv (ncbi-geneid, uniprot, etc.) --------
    for raw_id in gene_ids:
        clean_id = limpiar_id_gen(raw_id)
        source_id = f"{source_db}:{clean_id}"

        try:
            logging.info("Convirtiendo %s -> KEGG...", source_id)
            handle = REST.kegg_conv(org_code, source_id)
            txt = handle.read().strip()
            handle.close()
        except Exception as e:
            logging.warning("Error al consultar kegg_conv para %s: %s", source_id, e)
            time.sleep(delay)
            continue

        if not txt:
            logging.warning("Sin resultado de kegg_conv para %s", source_id)
            time.sleep(delay)
            continue

        kegg_id = ""
        for line in txt.splitlines():
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            src, dst = parts[0].strip(), parts[1].strip()
            if src == source_id:
                kegg_id = dst
                break

        if not kegg_id:
            logging.warning(
                "No se encontró un ID de KEGG para %s (organismo %s)",
                source_id,
                org_code,
            )
            time.sleep(delay)
            continue

        registros.append(
            {
                "gene_id": clean_id,
                "source_id": source_id,
                "kegg_id": kegg_id,
            }
        )

        time.sleep(delay)

    df_ids = pd.DataFrame(registros)
    logging.info("Conversión a KEGG completada para %d genes", df_ids.shape[0])
    return df_ids

# ---------------------------------------------------------------------
# Consulta y parseo de anotación con kegg_get (en chunks)
# ---------------------------------------------------------------------


def parsear_entrada_kegg(raw_entry: str) -> Dict[str, str]:
    """
    Parsea una entrada de texto plano devuelta por REST.kegg_get().
    Extrae: NAME, DEFINITION, KO, PATHWAY(s).
    """
    name = ""
    definition = ""
    ko = ""
    pathways: List[str] = []

    current_field = None

    for line in raw_entry.splitlines():
        if not line.strip():
            continue

        tag = line[:12].strip()
        data = line[12:].strip()

        if tag:
            current_field = tag
        else:
            tag = current_field

        if tag == "NAME":
            if name:
                name += " " + data
            else:
                name = data
        elif tag == "DEFINITION":
            if definition:
                definition += " " + data
            else:
                definition = data
        elif tag == "KO":
            ko = data.split()[0]
        elif tag == "PATHWAY":
            pathways.append(data)

    return {
        "name": name,
        "definition": definition,
        "ko": ko,
        "pathways": "; ".join(pathways) if pathways else "",
    }


from urllib.error import HTTPError

def anotar_genes_con_kegg(
    kegg_ids: List[str],
    delay: float = 4.0,
    chunk_size: int = 10,
) -> pd.DataFrame:
    """
    Consulta KEGG para una lista de IDs de KEGG y devuelve anotación.

    Se consultan varios IDs a la vez en "chunks" usando kegg_get.
    Si KEGG responde 404 para un chunk, se hace un fallback a
    consultas individuales por ID.

    Parameters
    ----------
    kegg_ids : list of str
        Lista de IDs de KEGG (ej. 'eco:b0001', 'spo:SPAC1002.07c', ...).
    delay : float
        Tiempo de espera (segundos) entre peticiones.
    chunk_size : int
        Número máximo de IDs por petición.

    Returns
    -------
    pandas.DataFrame
        Tabla con columnas:
        - kegg_id
        - name
        - definition
        - ko
        - pathways
    """
    registros: List[Dict[str, str]] = []

    # quitar duplicados conservando orden
    unique_ids = list(dict.fromkeys(kegg_ids))

    def procesar_raw(raw: str, org_prefix: str) -> None:
        """Parsea el texto devuelto por kegg_get y agrega a 'registros'."""
        if not raw or not raw.strip():
            return

        for entry in raw.strip().split("///"):
            entry = entry.strip()
            if not entry:
                continue

            # Recuperar el ID desde la línea ENTRY
            entry_id_raw = ""
            for line in entry.splitlines():
                tag = line[:12].strip()
                data = line[12:].strip()
                if tag == "ENTRY":
                    # Puede ser 'b0002' o 'spo:SPAC1002.07c'
                    entry_id_raw = data.split()[0]
                    break

            if not entry_id_raw:
                continue

            # Si no trae prefijo, se lo agregamos (eco:b0002, spo:SPAC1002.07c)
            if ":" in entry_id_raw:
                entry_id_full = entry_id_raw
            else:
                entry_id_full = f"{org_prefix}:{entry_id_raw}"

            info = parsear_entrada_kegg(entry)
            info["kegg_id"] = entry_id_full
            registros.append(info)

    for start in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[start:start + chunk_size]
        if not chunk:
            continue

        org_prefix = chunk[0].split(":", 1)[0]

        logging.info(
            "Consultando KEGG (kegg_get) para IDs %d–%d de %d...",
            start + 1,
            start + len(chunk),
            len(unique_ids),
        )

        # --- 1) Intento: pedir el chunk completo como lista ---
        try:
            handle = REST.kegg_get(chunk)  # Biopython acepta lista de IDs
            raw = handle.read()
            handle.close()
            if raw.strip():
                procesar_raw(raw, org_prefix)
                time.sleep(delay)
                continue
            else:
                logging.warning(
                    "KEGG no devolvió contenido para el chunk %d–%d",
                    start + 1,
                    start + len(chunk),
                )
        except HTTPError as e:
            logging.warning(
                "Error HTTP al consultar kegg_get para el chunk %d–%d: %s",
                start + 1,
                start + len(chunk),
                e,
            )
        except Exception as e:
            logging.warning(
                "Error al consultar kegg_get para el chunk %d–%d: %s",
                start + 1,
                start + len(chunk),
                e,
            )

        # --- 2) Fallback: intentar uno por uno en este chunk ---
        #logging.info(
            #"Intentando fallback (kegg_get por gen) para el chunk %d–%d...",
            #start + 1,
           # start + len(chunk),
       # )
        for kid in chunk:
            try:
                h2 = REST.kegg_get(kid)
                raw_single = h2.read()
                h2.close()
            except Exception as e2:
                logging.warning(
                    "Error al consultar kegg_get para %s: %s",
                    kid,
                    e2,
                )
                continue

            if not raw_single.strip():
                logging.warning("KEGG no devolvió contenido para %s", kid)
                continue

            procesar_raw(raw_single, org_prefix)
            time.sleep(delay)

    df_ann = pd.DataFrame(registros)
    if not df_ann.empty:
        df_ann = df_ann.drop_duplicates(subset="kegg_id", keep="first")

    logging.info("Anotación KEGG obtenida para %d IDs únicos", df_ann.shape[0])
    return df_ann


# ---------------------------------------------------------------------
# Pipeline: lista de genes -> conversión -> anotación
# ---------------------------------------------------------------------


def anotar_lista_genes_con_kegg(
    genes_file: str,
    org_code: str,
    source_db: str,
    set_label: str | None = None,
) -> pd.DataFrame:
    """
    Pipeline completo para una lista de genes ya filtrados (UP o DOWN):

    1. Lee la lista de genes.
    2. Convierte IDs a IDs de KEGG.
    3. Consulta KEGG para esos IDs en chunks.
    4. Devuelve un DataFrame combinado.
    """
    genes_raw = cargar_lista_genes(genes_file)
    if not genes_raw:
        raise ValueError(f"La lista de genes está vacía: {genes_file}")

    df_ids = convertir_ids_a_kegg(
        gene_ids=genes_raw,
        source_db=source_db,
        org_code=org_code,
    )

    if df_ids.empty:
        logging.warning(
            "No se obtuvo ningún ID de KEGG. Revisa 'source_db', 'org_code' "
            "y el formato de los IDs de la lista."
        )
        return df_ids

    unique_kegg_ids = df_ids["kegg_id"].unique().tolist()
    df_kegg = anotar_genes_con_kegg(unique_kegg_ids)

    df_final = df_ids.merge(df_kegg, on="kegg_id", how="left")

    if set_label is not None:
        df_final["set_label"] = set_label

    return df_final


# ---------------------------------------------------------------------
# Main con argparse
# ---------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    """
    Define y parsea los argumentos de línea de comando.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Anotar funcionalmente una LISTA de genes (UP/DOWN) usando KEGG. "
            "Recibe un archivo de texto con IDs de genes (uno por línea)."
        )
    )

    parser.add_argument(
        "--genes-file",
        required=True,
        help="Ruta al archivo de texto con la lista de genes (uno por línea).",
    )
    parser.add_argument(
        "--org-code",
        required=True,
        help="Código de organismo en KEGG (ej. eco, spo, phavu).",
    )
    parser.add_argument(
        "--source-db",
        required=True,
        help=(
            "Base de datos origen para kegg_conv (ej. ncbi-geneid, uniprot). "
            "Si los IDs de la lista ya son los usados por KEGG (locus_tag, "
            "por ejemplo b0001 para E. coli), usa 'direct'."
        ),
    )
    parser.add_argument(
        "--set-label",
        default=None,
        help=(
            "Etiqueta opcional para la lista de genes (ej. 'UP', 'DOWN'). "
            "Se agrega como columna 'set_label' en la salida."
        ),
    )
    parser.add_argument(
        "--out",
        default="kegg_annotation.tsv",
        help="Archivo de salida con la anotación KEGG (TSV).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Nivel de logging (por defecto INFO).",
    )

    return parser.parse_args()


def main() -> None:
    """
    Punto de entrada principal cuando se ejecuta como script.
    """
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
    )

    logging.info("Iniciando anotación KEGG para lista de genes...")

    df_final = anotar_lista_genes_con_kegg(
        genes_file=args.genes_file,
        org_code=args.org_code,
        source_db=args.source_db,
        set_label=args.set_label,
    )

    # Guardamos siempre en TSV
    df_final.to_csv(args.out, sep="\t", index=False)
    logging.info("Archivo de anotación KEGG guardado en: %s", args.out)


if __name__ == "__main__":
    main()
