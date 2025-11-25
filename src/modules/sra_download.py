"""
sra_download.py
Autor: Alondra Márquez
Fecha: 2025-10-28

Funciones auxiliares para trabajar con datos de SRA:

- Verificar que archivos FASTQ existan y no estén vacíos.
- Verificar que un archivo tiene formato FASTQ comprimido válido.
- Descargar SRR usando `prefetch` y convertirlos a FASTQ con `fastq-dump`,
  organizando cada SRR en su propio directorio.
- Concatenar múltiples archivos FASTQ (por muestra) en un único par de
  archivos FASTQ comprimidos (R1/R2).

Este módulo está pensado para ser utilizado por el script principal
`main_download.py` 

"""

import subprocess
import os
import logging
import gzip

#Funcion creada para verificar que el archivo exista y que no este vacio 
def check_file(file_path, min_size=1000):
    """
    Verifica que un archivo exista y no esté vacío (por debajo de un tamaño mínimo).

    Returns
    -------
    bool
        True si el archivo existe y tiene al menos `min_size` bytes;
        False en caso contrario.
    """
    if not os.path.exists(file_path):
        logging.error(f"Archivo no encontrado: {file_path}")
        return False
    if os.path.getsize(file_path) < min_size:
        logging.warning(f"Archivo muy pequeño o incompleto: {file_path}")
        return False
    return True


#Comprueba que el archivo FASTQ esté comprimido y extension correcta
def verify_fastq(file_path):
    """
    Comprueba que un archivo sea un FASTQ gzip válido.

    La verificación incluye:
    - Que la extensión termine en ".fastq.gz".
    - Que el archivo exista y no esté vacío (`check_file`).
    - Que la primera línea (header) comience con '@', patrón típico de FASTQ.

    Parameters
    ----------
    file_path : str
        Ruta al archivo FASTQ comprimido (.fastq.gz).

    Returns
    -------
    bool
        True si el archivo parece ser un FASTQ válido, False en caso contrario.
    """

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
def prefetch_fastqdump (srr_id, fastq_dir, retries=2): #retries: numero de reintentos si falla la descarga (default=2).
    """
    Descarga un SRR y lo convierte a FASTQ (single o paired) usando un script Bash.

    El flujo es:
    - Crear un directorio específico para el SRR dentro de `fastq_dir`.
    - Generar un script Bash con los comandos:
        * `prefetch` para descargar el SRR.
        * `fastq-dump --split-files --gzip` para generar los FASTQ.
    - Ejecutar el script, con un número de reintentos configurable.
    - Verificar qué archivos FASTQ existen (paired o single).

    Parameters
    ----------
    srr_id : str
        Identificador SRR de SRA (por ejemplo, "SRR10192862").
    fastq_dir : str
        Directorio base donde se crearán las carpetas para cada SRR.
    retries : int, optional
        Número máximo de reintentos si la descarga falla (por defecto: 2).

    Returns
    -------
    (str or None, str or None)
        Tupla con rutas a los archivos FASTQ generados:
        - (file_r1, file_r2) si los datos son paired-end.
        - (file_single, None) si los datos son single-end.
        - (None, None) si no se pudo obtener un FASTQ válido.
    """
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
    elif verify_fastq(file_s):  # single-end clásico
        logging.info(f"{srr_id}: Lecturas single-end (archivo único)")
        return file_s, None
    elif verify_fastq(file_r1):
        logging.info(f"{srr_id}: Lecturas single-end ")
        return file_r1, None
    else:
        logging.error(f"No se encontraron FASTQ válidos para {srr_id}")
        return None, None


#Concatenar archivos FASTQ y comprimir 
def concatenate(fastq_dir, files_r1, files_r2=None, sample_name="sample"):
     
    """
    Concatena múltiples archivos FASTQ (R1 y opcionalmente R2) en archivos
    finales comprimidos por muestra.

    El resultado se guarda en:
        <fastq_dir>/<sample_name>/<sample_name>_1.fastq.gz
        <fastq_dir>/<sample_name>/<sample_name>_2.fastq.gz  (si files_r2 no es None)

    Parameters
    ----------
    fastq_dir : str
        Directorio base donde se creará el directorio de la muestra.
    files_r1 : list of str
        Lista de rutas a los archivos FASTQ de lecturas forward (R1).
    files_r2 : list of str or None, optional
        Lista de rutas a los archivos FASTQ de lecturas reverse (R2),
        o None si los datos son single-end.
    sample_name : str, optional
        Nombre de la muestra; se utilizará para nombrar el directorio y
        los archivos de salida (por defecto: "sample").

    Raises
    ------
    FileNotFoundError
        Si alguno de los archivos de entrada no existe.
    RuntimeError
        Si ocurre algún error durante la concatenación con `zcat`/`gzip`.
    """
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

