'''
main_download.py
Autor: Alondra Marquez
Fecha: 2025-10-28
Descripción: 

'''

from Bio import Entrez
import pandas as pd
from io import StringIO
import logging


# Buscar el UID del BioProject
#Busca el UID de un BioProject usando su accession. Si hay error , se registra y retorna None.
def set_email(email):
    Entrez.email = email
def consulta_bioproject(project_accession):
    
    try:
        logging.info(f"Iniciando consulta BioProject: {project_accession}")
        handle = Entrez.esearch(db="bioproject", term=project_accession)
        search_results = Entrez.read(handle)
        handle.close()
        
        if not search_results.get("IdList"):
            logging.error(f"No se encontró el BioProject {project_accession}")
            return None
        
        project_uid = search_results["IdList"][0]
        logging.info(f"BioProject encontrado: UID {project_uid}")
        return project_uid
    
    except Exception as e:
        logging.error(f"Error consultando BioProject {project_accession}: {e}")
        return None


#Resumen del proyecto: con el IUD y con la funcion esummary de Biopython sobre la base de datos bioproject se puede obtener el resumen del proyecto 
def bioproject_summary(project_uid):
    try:
        logging.info(f"Obteniendo resumen del BioProject UID {project_uid}")
        handle = Entrez.esummary(db="bioproject", id=project_uid)
        project_record = Entrez.read(handle)
        handle.close()
        
        doc_summary = project_record.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        if not doc_summary:
            logging.error(f"No se encontró información para UID {project_uid}")
            return None
    
        project = doc_summary[0]
        summary = {
            "Acceso": project.get("Project_Accession"),
            "Título": project.get("Project_Title"),
            "Organismo": project.get("Organism_Name"),
            "Descripción": project.get("Project_Description", "No disponible"),
        }
        logging.info(f"Resumen obtenido para UID {project_uid}")
        return summary
    except Exception as e:
        logging.error(f"Error obteniendo resumen de BioProject UID {project_uid}: {e}")
        return None


# Relacionar BioProject con GEO (GSE) con el uid del bioproject
def links_bioproject(project_uid, project_accession):
    try:
        logging.info(f"Obteniendo enlaces GEO para BioProject UID {project_uid}")
        handle = Entrez.elink(dbfrom="bioproject", db="gds", id=project_uid)
        linksBioProj = Entrez.read(handle)
        handle.close()
        
        linksets = linksBioProj[0].get("LinkSetDb", [])
        if linksets:
            geo_ids = [link["Id"] for link in linksets[0].get("Link", [])]
            logging.info(f"{len(geo_ids)} series GEO asociadas a {project_accession}")
            return geo_ids
        else:
            logging.info(f"No se encontraron enlaces GEO para {project_accession}")
            return []
    except Exception as e:
        logging.error(f"Error obteniendo enlaces GEO para {project_accession}: {e}")
        return []


def geo_summary(geo_uid):
    # Verificando la información de GEO
    try:
        logging.info(f"Obteniendo resumen de GEO UID {geo_uid}")
        handle = Entrez.esummary(db="gds", id=geo_uid)
        summary = Entrez.read(handle)
        handle.close()
        
        if not summary:
            logging.error(f"No se encontró información para GEO UID {geo_uid}")
            return None
        
        geo_record = summary[0]
        data = {
            "Accession": geo_record.get("Accession"),
            "Título": geo_record.get("title"),
            "Resumen": geo_record.get("summary"),
            "Organismo": geo_record.get("taxon"),
            "N° muestras": geo_record.get("n_samples"),
            "PubMed IDs": geo_record.get("PubMedIds"),
            "Samples": []
        }

        samples_list = geo_record.get("Samples", [])
        for sample in samples_list:
            data["Samples"].append({
                "GSM": sample.get("Accession"),
                "Título": sample.get("Title")
            })

        logging.info(f"Resumen de GEO UID {geo_uid} obtenido con {len(data['Samples'])} muestras")
        return data
    except Exception as e:
        logging.error(f"Error obteniendo resumen de GEO UID {geo_uid}: {e}")
        return None
  
# Devuelve los SRR asociados a una muestra GSM
def srr_(gsm_id):
    #Busca los SRR asociados a una muestra GSM en el SRA
    try:
        logging.info(f"Buscando SRR para GSM {gsm_id}")
        #Busca los SRR asociados a una muestra GSM en el SRA
        handle = Entrez.esearch(db="sra", term=gsm_id)
        record = Entrez.read(handle)
        handle.close()
        
        if not record.get("IdList"):
            logging.warning(f"No se encontraron SRR para {gsm_id}")
            return pd.DataFrame()
        
        sra_ids = record["IdList"]
        handle = Entrez.efetch(db="sra", id=sra_ids, rettype="runinfo", retmode="text")
        runinfo = handle.read()
        handle.close()
        
        if isinstance(runinfo, bytes):
            runinfo = runinfo.decode("utf-8")

        df = pd.read_csv(StringIO(runinfo))
        required_cols = ["Run", "Experiment", "Platform", "LibraryName", "LibraryLayout",
                         "Sample", "ScientificName", "SampleName"]
        
        return df[required_cols]   
    except Exception as e:
        logging.error(f"Error consultando SRR para GSM {gsm_id}: {e}")
        return pd.DataFrame()
