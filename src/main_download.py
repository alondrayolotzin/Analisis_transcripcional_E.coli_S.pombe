"""
main_download.py
Autor: Alondra Márquez
Fecha: 2025-10-28

Descripción
-----------
Pipeline para:

1. Consultar un BioProject en NCBI y obtener información básica.
2. Recuperar el GEO asociado y listar las muestras (GSM) y sus SRR.
3. Descargar los SRR (lecturas crudas), organizarlos por especie y por muestra
   (GSM) y concatenar los FASTQ por muestra.
4. Generar una tabla de metadata con información mínima por muestra
   (sample, species, path).
5. Descargar los genomas de referencia y anotaciones (fna/gff/faa) para
   los organismos de interés.

Uso de ejemplo
--------------
python main_download.py \
    --bioproject PRJNA574477 \
    --fastq-dir /export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/data \
    --genome-dir /export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/data \
    --organismos "Schizosaccharomyces pombe" "Escherichia coli" \
    --accessions GCF_000002945.2 GCF_000005845.2 \
    --incluir fna gff \
    --log /export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/logs/pipeline_master.log \
    --email alondram@lcg.unam.mx
"""


import argparse, logging, os, csv 
from modules import bioproject, sra_download, genome_reference

def main():

    #  # === Definición y parseo de argumentos de línea de comandos ===
    parser = argparse.ArgumentParser(description="Pipeline completo de consulta, descarga de datos SRR desde BioProject y descarga de genomas de referencia.")
    parser.add_argument("--bioproject", required=True, help="Accession del BioProject (e.g., PRJNA574477)")
    parser.add_argument("--fastq-dir", required=True, help="Directorio de salida donde se guardarán los FASTQ")
    parser.add_argument("--genome-dir", required=True, help="Directorio de salida donde se guardarán los genomas de referencia") # nargs Indica que el argumento acepta una o más palabras o valores (separados por espacios)
    parser.add_argument("--organismos", nargs="+", required=True, help="Lista de organismos para descargar genomas") 
    parser.add_argument("--accessions", nargs="+", required=True, help="Lista de accesiones correspondientes a los organismos")
    parser.add_argument("--incluir", nargs="+", default=["fna", "gff"], choices=["fna", "gff", "faa"], help="Tipos de archivos a descargar")
    parser.add_argument("--retmax", type=int, default=50, help="Número máximo de ensamblados a buscar (por defecto: 50).")
    parser.add_argument("--log", default="pipeline_master.log", help="Ruta al archivo log general")
    parser.add_argument("--email", required=True, help="Correo electrónico para consultas Entrez")
    args = parser.parse_args()

    if len(args.organismos) != len(args.accessions):
        raise ValueError("El numero de organismos y el numero de accesiones debe coincidir.")
    
   
    # === Configuración de Entrez para los módulos que lo requieren ===
    bioproject.set_email(args.email)
    genome_reference.set_email(args.email)

     # === Preparación de directorios de salida ===
    fastq_dir = args.fastq_dir
    genome_dir = args.genome_dir

    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(genome_dir, exist_ok=True)
    log_dir = os.path.dirname(args.log)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

     # === Configuración de logging (archivo + consola) ===
    log_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Handler para archivo
    file_handler = logging.FileHandler(args.log)
    file_handler.setFormatter(log_formatter)

    # Handler para consola
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)

    # Limpiar handlers previos y registrar los nuevos
    logging.getLogger().handlers = []  # limpia cualquier handler previo
    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(file_handler)
    logging.getLogger().addHandler(console_handler)

    # Mensaje de inicio del pipeline
    logging.info("=== Inicio del pipeline de BioProject ===")

    
    # === Consulta del BioProject ===
    logging.info("\n=== Consulta BioProject ===")
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
    print(f"Descripción: {information.get('Descripción', 'No disponible')}\n")

    
    # === Recuperar GEO asociado al BioProject ===

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
    
     # === Consulta de SRR por muestra (GSM) ===

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
        
        all_samples_srr[gsm_id] = df_srr  #df_srr["Run"].tolist()

    if not all_samples_srr:
        logging.error("No se encontraron SRR válidos para ninguna muestra. Finalizacion de la ejecucion" )
        return

    # === Descarga de SRR, concatenación de FASTQ y construcción de metadata ===
    # Cada fila de metadata: [sample (GSM), species, path_a_R1]
    metadata = []  # filas: [sample, species, path]

    for gsm_id, df_srr in all_samples_srr.items():

        # Determinar especie a partir de ScientificName 
        species = df_srr["ScientificName"].iloc[0].replace(" ", "_")
        species_dir = os.path.join(fastq_dir, species)
        os.makedirs(species_dir, exist_ok=True)

        r1_files = []
        r2_files = []

        # Descargar cada SRR asociado a esa muestra
        for srr_id in df_srr["Run"].tolist():
            logging.info(f"Iniciando descarga de {srr_id} para {gsm_id}")
            r1, r2 = sra_download.prefetch_fastqdump(srr_id, species_dir)

            if r1: 
                r1_files.append(r1) 
            if r2: 
                r2_files.append(r2)

        # Concatenar SRR por muestra y registrar en metadata
        if r1_files:
            try:
                sra_download.concatenate(
                    species_dir,
                    r1_files, 
                    r2_files if r2_files else None, 
                    sample_name=gsm_id
                )
                fastq_path = os.path.join(species_dir, gsm_id, f"{gsm_id}_1.fastq.gz")
                metadata.append([gsm_id, species, fastq_path])
            except Exception as e:
                logging.error(f"Error al concatenar FASTQ de {gsm_id}: {e}")
        else:
            logging.warning(f"No se descargaron archivos válidos para {gsm_id}. Se omitira concatenacion.")
        
    # === Guardar metadata para análisis posteriores (e.g. DESeq2) ===
    metadata_dir = os.path.join(fastq_dir, "metadata")
    os.makedirs(metadata_dir, exist_ok=True) 

    # Ruta final del CSV dentro de data/metadata/
    metadata_file = os.path.join(metadata_dir, "samples.csv")

    with open(metadata_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "species", "path"])
        writer.writerows(metadata)

    logging.info(f"Metadata generada en {metadata_file}")
    

     # === Descarga de genomas de referencia para cada organismo ===
        #zip() empareja elemento por elemento entre dos listas.
    for org, acc in zip(args.organismos, args.accessions): #itera sobre las dos listas generadas para le nombre de los organimos y para las accesiones 
        # Carpeta específica por organismo + accession, para no mezclar genomas
        subdir = os.path.join(args.genome_dir, f"{org.replace(' ', '_')}_{acc}") #ruta donde se guardaran los archivos de cada organismo 
        os.makedirs(subdir, exist_ok=True) #crea los subdirectorios 
        logging.info(f"Buscando ensambles de {org} (accession solicitada: {acc})")

        try:
            df = genome_reference.assemblies_por_organismo(org, retmax=args.retmax)#funcion de busqueda de ensamblados por organismo en NCBI, devuelve data frame 
        except Exception as e:
            logging.error(f"Fallo al buscar ensamblados para {org}: {e}")
            continue

        try:
            genome_reference.descargar_ensamblado(df, acc, subdir, incluir=args.incluir) #usa el DataFrame para localizar la fila con la accessión acc y descargar los archivos (fna/gff/faa según args.incluir) al subdir 
            logging.info(f"Descarga completada para {org}: {acc}")
        except ValueError as e:
            logging.error(f"Accessión {acc} no encontrada para {org}: {e}")

    logging.info("=== Pipeline completado ===")
    print("\nPipeline completado. Revisar logs para detalles.")

      
if __name__ == "__main__":
    main()

