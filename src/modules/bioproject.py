from Bio import Entrez
import pandas as pd
from io import StringIO
import os

#
Entrez.email = "alondram@chaac.lcg.unam.mx"

def consulta_bioproject(project_accesion):
# Buscar el UID del BioProject
    handle = Entrez.esearch(db="bioproject", term=project_accesion)
    search_results = Entrez.read(handle)
    handle.close()

    #Regresar el IUD del bioproject 
    if not search_results["IdList"]:
        raise ValueError(f"No se encontro el BioProject {project_accesion}")

    project_uid= search_results["IdList"][0]
    return project_uid

#Resumen del proyecto: con el IUD y con la funcion esummary de Biopython sobre la base de datos bioproject se puede obtener el resumen del proyecto 
def bioproject_summary(project_uid):
    handle = Entrez.esummary(db="bioproject", id=project_uid)
    project_record = Entrez.read(handle)
    handle.close()

    # Seleccionando la información relacionada al proyecto
    project = project_record["DocumentSummarySet"]["DocumentSummary"][0]
 
    # Accediendo a campos especificos del diccionario
    summary = {
        "Acceso": project.get("Project_Accession"),
        "Título": project.get("Project_Title"),
        "Organismo": project.get("Organism_Name"),
        "Descripción:": project.get("Project_Description", "No disponible"),
    }
    return summary

# Relacionar BioProject con GEO (GSE) con el uid del bioproject
def links_bioproject(project_uid, project_accesion):

    handle = Entrez.elink(dbfrom="bioproject", db="gds", id=project_uid)
    linksBioProj = Entrez.read(handle)
    handle.close()
    
    if linksBioProj[0]["LinkSetDb"]:
        geo_ids = [link["Id"] for link in linksBioProj[0]["LinkSetDb"][0]["Link"]]
        print(f"El proyecto {project_accesion} está asociado a GEO ({len(geo_ids)} series).")
        return geo_ids
    else:
        print(f"El proyecto {project_accesion} no tiene enlaces directos a GEO.")
        return []


def geo_summary(geo_uid):
    # Verificando la información de GEO
    handle = Entrez.esummary(db="gds", id=geo_uid)
    summary = Entrez.read(handle)
    handle.close()
   
    geo_record=summary[0]
    data = {
        "Accession": geo_record.get("Accession"),
        "Título": geo_record.get("title"),
        "Resumen": geo_record.get("summary"),
        "Organismo": geo_record.get("taxon"),
        "N° muestras": geo_record.get("n_samples"),
        "PubMed IDs": geo_record.get("PubMedIds"),
        "Samples": []
    }

    if "Samples" in geo_record:
        for sample in geo_record["Samples"]:
            data["Samples"].append({
                "GSM": sample["Accession"],
                "Título": sample["Title"]
            })

    return data
  
# Devuelve los SRR asociados a una muestra GSM
def srr_(gsm_id):
    #Busca los SRR asociados a una muestra GSM en el SRA
    handle = Entrez.esearch(db="sra", term=gsm_id)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        print(f"No se encontraron SRR para {gsm_id}")
        return pd.DataFrame()
    else:
        sra_ids = record["IdList"]

    handle = Entrez.efetch(db="sra", id=sra_ids, rettype="runinfo", retmode="text")
    runinfo = handle.read()
    handle.close()

    if isinstance(runinfo, bytes):
        runinfo = runinfo.decode("utf-8")

    df = pd.read_csv(StringIO(runinfo))
    return df[["Run", "Experiment", "Platform", "LibraryName", "LibraryLayout",
               "Sample", "ScientificName", "SampleName"]]
