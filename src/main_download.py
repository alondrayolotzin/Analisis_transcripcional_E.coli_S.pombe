'''
main_download.py
Autor: Alondra Marquez
Fecha: 2025-10-28
Descripción: 

Uso: 
python main_download.py \
    --bioproject PRJNA574477 \
    --fastq-dir /export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/data \
    --genome-dir /export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/data \
    --organismos "Schizosaccharomyces pombe" "Escherichia coli" \
    --accessions GCF_000002945.2 GCF_000005845.2 \
    --incluir fna gff \
    --email alondram@lcg.unam.mx
PRJNA574477
/export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/data


'''


import argparse, logging, os
from modules import bioproject, sra_download, genome_reference

def main():

    # === Argumentos del pipeline ===
    parser = argparse.ArgumentParser(description="Pipeline completo de consulta, descarga de datos SRR desde BioProject y descarga de genomas de referencia.")
    parser.add_argument("--bioproject", required=True, help="Accession del BioProject (e.g., PRJNA574477)")
    parser.add_argument("--fastq-dir", required=True, help="Directorio de salida donde se guardarán los FASTQ")
    parser.add_argument("--genome-dir", required=True, help="Directorio de salida donde se guardarán los genomas de referencia") # nargs Indica que el argumento acepta una o más palabras o valores (separados por espacios)
    parser.add_argument("--organismos", nargs="+", required=True, help="Lista de organismos para descargar genomas") 
    parser.add_argument("--accessions", nargs="+", required=True, help="Lista de accesiones correspondientes a los organismos")
    parser.add_argument("--incluir", nargs="+", default=["fna", "gff"], choices=["fna", "gff", "faa"],
                        help="Tipos de archivos a descargar")
    parser.add_argument("--retmax", type=int, default=50, help="Número máximo de ensamblados a buscar (por defecto: 50).")
    parser.add_argument("--log", default="pipeline_master.log", help="Ruta al archivo log general")
    parser.add_argument("--email", required=True, help="Correo electrónico para consultas Entrez")
    args = parser.parse_args()

    # Logging global del pipeline
    # --------------------------
   # Setea el email en los módulos que lo necesitan
    bioproject.set_email(args.email)
    genome_reference.set_email(args.email)

# Directorios de salida
    fastq_dir = args.fastq_dir
    genome_dir = args.genome_dir

    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(genome_dir, exist_ok=True)
    log_dir = os.path.dirname(args.log)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    # Configuración de logging dual (archivo + consola)
    log_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Handler para archivo
    file_handler = logging.FileHandler(args.log)
    file_handler.setFormatter(log_formatter)

    # Handler para consola
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)

    # Asigna ambos handlers al root logger
    logging.getLogger().handlers = []  # limpia cualquier handler previo
    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(file_handler)
    logging.getLogger().addHandler(console_handler)

    # Mensaje de inicio del pipeline
    logging.info("=== Inicio del pipeline de BioProject ===")

    print("\n=== Consulta BioProject ===")
    project_uid = bioproject.consulta_bioproject(args.bioproject)
    if not project_uid:
        logging.error(f"BioProject {args.bioproject} no encontrado. Saliendo del pipeline.")
        return
    information = bioproject.bioproject_summary(project_uid)
    if not information:
        logging.error(f"No se pudo obtener resumen de {args.bioproject}. Saliendo del pipeline.")
        return

    print("\n=== Información del BioProject ===")
    print(f"Acceso: {information.get('Acceso', 'No disponible')}")
    print(f"Título: {information.get('Título', 'No disponible')}")
    print(f"Organismo: {information.get('Organismo', 'No disponible')}")
    print(f"Descripción: {information.get('Descripción:', 'No disponible')}\n")

    # Obtener GEO IDs

    geo_ids = bioproject.links_bioproject(project_uid, args.bioproject)
    if not geo_ids:
        logging.warning(f"No hay GEO asociado a {args.bioproject}. Finalizando pipeline.")
        return
    
    geo_info = bioproject.geo_summary(geo_ids[0])
    if not geo_info:
        logging.error(f"No se pudo obtener información de GEO UID {geo_ids[0]}.")
        return

    print("\n=== Información GEO ===")
    print(f"Accession: {geo_info.get('Accession', 'No disponible')}")
    print(f"Título: {geo_info.get('Título', 'No disponible')}")
    print(f"Resumen: {geo_info.get('Resumen', 'No disponible')}")
    print(f"Organismo: {geo_info.get('Organismo', 'No disponible')}")
    print(f"Número de muestras: {geo_info.get('N° muestras', 'No disponible')}")
    print(f"PubMed IDs: {geo_info.get('PubMed IDs', 'No disponible')}\n")

    print("=== Muestras asociadas ===")
    for sample in geo_info.get("Samples", []):
        print(f"{sample['GSM']}: {sample['Título']}")
    
    # SRR de cada muestra
    print("\n=== SRR asociados a cada muestra ===")
    all_samples_srr = {} 

    for sample in geo_info["Samples"]:
        gsm_id = sample["GSM"]
        df_srr = bioproject.srr_(gsm_id)

        if df_srr.empty:
            logging.warning(f"No se encontraron SRR para {gsm_id}.")
            continue

        print(f"\nMuestra {gsm_id} - {sample['Título']}")
        print(df_srr[["Run", "Experiment", "Platform", "LibraryName", 
                      "LibraryLayout", "Sample", "ScientificName", "SampleName"]].to_string(index=False))
        
        all_samples_srr[gsm_id] = df_srr["Run"].tolist()

     # Descarga y concatenación SRR

    for gsm_id, srr_list in all_samples_srr.items():
        r1_files = []
        r2_files = []

        for srr_id in srr_list:
            logging.info(f"Iniciando descarga de {srr_id} para {gsm_id}")
            r1, r2 = sra_download.prefetch_fastqdump(srr_id, args.fastq_dir)

            if r1: 
                r1_files.append(r1) 
            if r2: 
                r2_files.append(r2)

        if r1_files:
            sra_download.concatenate(args.fastq_dir, r1_files, r2_files if r2_files else None, sample_name=gsm_id)
        else:
            logging.warning(f"No se descargaron archivos válidos para {gsm_id}. Se omitira concatenacion.")

      # === Descarga de genomas de referencia ===
        #zip() empareja elemento por elemento entre dos listas.
    for org, acc in zip(args.organismos, args.accessions): #itera sobre las dos listas generadas para le nombre de los organimos y para las accesiones 
        subdir = os.path.join(args.genome_dir, org.replace(" ", "_")) #ruta donde se guardaran los archivos de cada organismo 
        os.makedirs(subdir, exist_ok=True) #crea los subdirectorios 
        logging.info(f"Buscando ensambles de {org} (accession solicitada: {acc})")
        try:
            df = genome_reference.assemblies_por_organismo(org, retmax=args.retmax)#funcion de busqueda de ensamblados por organismo en NCBI, devuelve data frame 
        except Exception as e:
            logging.error(f"Fallo al buscar ensamblados para {org}: {e}")
        continue

        try:
            genome_reference.descargar_ensamblado(df, acc, subdir, incluir=args.incluir) #usa el DataFrame para localizar la fila con la accessión acc y descargar los archivos (fna/gff/faa según args.incluir) al subdir 
            logging.info(f"Descarga completada para {org}: {names}")
        except ValueError as e:
            logging.error(f"Accessión {acc} no encontrada para {org}: {e}")

    logging.info("=== Pipeline completado ===")
    print("\nPipeline completado. Revisar logs para detalles.")

      
if __name__ == "__main__":
    main()

