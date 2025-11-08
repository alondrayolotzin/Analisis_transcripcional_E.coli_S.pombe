'''
main_download.py
Autor: Alondra Marquez
Fecha: 2025-10-28
Descripción: 


'''

import subprocess
import os
import logging
import gzip

#Funcion creada para verificar que el archivo exista y que no este vacio 
def check_file(file_path, min_size=1000):
    if not os.path.exists(file_path):
        logging.error(f"Archivo no encontrado: {file_path}")
        return False
    if os.path.getsize(file_path) < min_size:
        logging.warning(f"Archivo muy pequeño o incompleto: {file_path}")
        return False
    return True


#Comprueba que el archivo FASTQ esté comprimido y extension correcta
def verify_fastq(file_path):
    if not file_path.endswith((".fastq.gz")):
        logging.warning(f"{file_path} no tiene una extensión FASTQ válida.")
        return False
    if not check_file(file_path):
        return False
    try:
        with gzip.open(file_path, "rb") as f:
            header = f.readline()
            if not header.startswith(b"@"): # se usa @ porque es un patron tipico de los headers de este tipo de archivo
                logging.warning(f"{file_path} no parece ser un archivo FASTQ válido (sin '@').")
                return False
        return True
    except Exception as e:
        logging.error(f"Error al verificar {file_path}: {e}")
        return False
    
#Usa Bash para descargar y convertir un SRR en FASTQ dentro de su propio directorio.
def  prefetch_fastqdump (srr_id, fastq_dir, retries=2): #retries: numero de reintentos si falla la descarga (default=2).

    #rutas para almacenar las descagas de los srr dentro del directorio fastq_dir y su log especifico 
    srr_dir = os.path.join(fastq_dir, srr_id)
    os.makedirs(srr_dir, exist_ok=True)

    bash_script = os.path.join(srr_dir, f"download_{srr_id}.sh") #Ruta donde se va a escribir el script bash que contiene los comandos de descarga/conversion    
    with open(bash_script, "w") as f:
        f.write("#!/bin/bash\nset -e\n")
        f.write(f"echo '=== Descargando {srr_id} ==='\n") #(#!/bin/bash) y set -e para detenerse si hay errores.
        f.write(f"prefetch --output-directory {srr_dir} {srr_id}\n") #prefetch para descargar el SRR a srr_dir
        f.write(f"fastq-dump --split-files --gzip -O {srr_dir} {srr_id}\n") #fastq-dump para convertir a FASTQ, separando pares (--split-files) y comprimiendo (--gzip).

    os.chmod(bash_script, 0o755) #cambio permisos del script para que sea ejecutable 

    logging.info(f"Script Bash generado: {bash_script}")#Registra en el log que el script fue creado para el srr

    for intento in range(1, retries + 1):
        try:
            subprocess.run(["bash", bash_script], check=True, capture_output=True, text=True) #subprocess.run ejecuta el script
            logging.info(f" Descarga completada del {srr_id}")
            break
        except subprocess.CalledProcessError as e:
            logging.error(f"Intento {intento} fallido para {srr_id}: {e.stderr}")
            if intento == retries:
                logging.error(f" Fallo definitivo de la descarga de {srr_id}")
                return None, None

    file_r1 = os.path.join(srr_dir, f"{srr_id}_1.fastq.gz") #Ruta esperada para R1
    file_r2 = os.path.join(srr_dir, f"{srr_id}_2.fastq.gz") #Ruta esperada para R2
    file_s = os.path.join(srr_dir, f"{srr_id}.fastq.gz") #Ruta esperada para single-end 

    if verify_fastq(file_r1) and verify_fastq(file_r2):
        logging.info(f"{srr_id}: Lecturas pair-end")
        return file_r1, file_r2
    elif verify_fastq(file_r1):
        logging.info(f"{srr_id}: Lecturas single-end ")
        return file_s, None
    else:
        logging.error(f"No se encontraron FASTQ válidos para {srr_id}")
        return None, None


#Concatenar archivos FASTQ y comprimir 
def concatenate(fastq_dir, files_r1, files_r2=None, sample_name="sample"):

    sample_dir = os.path.join(fastq_dir, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    logging.info(f"Iniciando concatenación para {sample_name}")

    #Recorre todos los archivos R1 y verifica que existan, si falta alguno, registra un error y lanza FileNotFoundError para detener el proceso
    for f in files_r1:
        if not os.path.exists(f):
            logging.error(f"Archivo faltante: {f}")
            raise FileNotFoundError(f"Archivo faltante: {f}")

    if files_r2:
        for f in files_r2:
            if not os.path.exists(f):
                logging.error(f"Archivo faltante: {f}")
                raise FileNotFoundError(f"Archivo faltante: {f}")
            
    #Define los nombres de salida finales para las lecturas forward y reverse, va a ser el nombre de la muestra 
    out_r1 = os.path.join(sample_dir, f"{sample_name}_1.fastq.gz")
    out_r2 = os.path.join(sample_dir, f"{sample_name}_2.fastq.gz") if files_r2 else None

    #comandos de bash para concatenar los archivos zcat descomprime los archivos en streamin, gzip vuelve a comprimir y se redirige salida a out_r1/out_r2.
    concatenar_r1 = f"zcat {' '.join(files_r1)} | gzip > {out_r1}"
    concatenar_r2 = f"zcat {' '.join(files_r2)} | gzip > {out_r2}" if files_r2 else None

    try:
        logging.info(f"Concatenando archivos R1 para {sample_name}")
        subprocess.run(concatenar_r1, shell=True, check=True) # permite pasar la cadena como comando de bash completo.
        logging.info(f"Archivo R1 concatenado en: {out_r1}") #regsitro en el log que la concatenacion fue exitosa 

        if files_r2:
            logging.info(f"Concatenando archivos R2 para {sample_name}")
            subprocess.run(concatenar_r2, shell=True, check=True)
            logging.info(f"Archivo R2 concatenado en: {out_r2}")

    except subprocess.CalledProcessError as e:
        logging.error(f"Error durante la concatenación: {e}")
        raise RuntimeError("Error en la concatenación de archivos FASTQ") from e
