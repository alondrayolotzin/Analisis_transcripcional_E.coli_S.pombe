from pathlib import Path
from Bio import Entrez
from modules.descarga import assemblies_por_organismo, descargar_ecoli_ref, descargar_spombe_ref
import argparse

def main(step: str):
    Entrez.email = "miryamzj@lcg.unam.mx"
    ROOT = Path(__file__).resolve().parent
    DATA = ROOT / "data"; DATA.mkdir(exist_ok=True)

    S_POMBE = "Schizosaccharomyces pombe"
    E_COLI  = "Escherichia coli"

    # Siempre podemos calcular estas tablas (son rápidas)
    df_pombe = assemblies_por_organismo(S_POMBE, 20)
    df_coli  = assemblies_por_organismo(E_COLI, 50)

    if step in ("assemblies", "all"):
        print("Assemblies S. pombe:", len(df_pombe))
        print("Assemblies E. coli :", len(df_coli))

    if step in ("download", "all"):
        out_p = DATA / "S_pombe_ASM294v3"
        out_c = DATA / "E_coli_MG1655"
        descargar_spombe_ref(df_pombe, outdir=out_p)   # cada función debería saltar si ya existe genome.fa
        descargar_ecoli_ref(df_coli,  outdir=out_c)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--step", choices=["assemblies","download","all"], default="all")
    args = p.parse_args()
    main(args.step)
