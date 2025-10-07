from modules import bioproject, sra_download

# Consulta de proyecto: se llama a el modulo biproject, que accede a Entrez
def main():
    project_accesion = "PRJNA574477"
    output_dir = "/export/storage/users/alondram/GENOMICAS/sem_2026/Biopython_sem26/Project_ecoli_spombe/data"

    print("\n=== Consulta BioProject ===")
    project_uid = bioproject.consulta_bioproject(project_accesion)
    information = bioproject.bioproject_summary(project_uid)


    print("\n=== Información del BioProject ===")
    print(f"Acceso: {information.get('Acceso', 'No disponible')}")
    print(f"Título: {information.get('Título', 'No disponible')}")
    print(f"Organismo: {information.get('Organismo', 'No disponible')}")
    print(f"Descripción: {information.get('Descripción:', 'No disponible')}\n")

    geo_ids = bioproject.links_bioproject(project_uid, project_accesion)
    if not geo_ids:
        return
    else:

        geo_info = bioproject.geo_summary(geo_ids[0])

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
            print(f"No se encontraron SRR para {gsm_id}")
            continue

        print(f"\nMuestra {gsm_id} - {sample['Título']}")

        print(df_srr[["Run", "Experiment", "Platform", "LibraryName", 
                      "LibraryLayout", "Sample", "ScientificName", "SampleName"]].to_string(index=False))
        
        all_samples_srr[gsm_id] = df_srr["Run"].tolist()

    for gsm_id, srr_list in all_samples_srr.items():
        srr_files_r1 = []
        srr_files_r2 = []

        for srr_id in srr_list:
            sra_download.prefetch_sra(srr_id, output_dir)
            r1, r2 = sra_download.fastq_dump(srr_id, output_dir)
            if r1: 
                srr_files_r1.append(r1) 
            if r2: 
                srr_files_r2.append(r2)
    
        sra_download.concatenate(output_dir, srr_files_r1, srr_files_r2, gsm_id)
        
if __name__ == "__main__":
    main()
