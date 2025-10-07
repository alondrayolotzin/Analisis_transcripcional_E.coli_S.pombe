import subprocess
import os

#Descarga los datos SRA usando prefetch
def prefetch_sra(srr_id, output_dir):
    print(f"Descarga de {srr_id} con prefetch")
    subprocess.run(["prefetch", "--output-directory", output_dir, srr_id], check=True)

#conversion de sra en FASTQ comprimidos
def fastq_dump(srr_id, output_dir): 
    subprocess.run(["fastq-dump", srr_id, "-O", output_dir,"--split-files", "--gzip"], check=True)

    # Comprobar si son paired-end o single-end
    file_r1 = os.path.join(output_dir, f"{srr_id}_1.fastq.gz")
    file_r2 = os.path.join(output_dir, f"{srr_id}_2.fastq.gz")
    file_s = os.path.join(output_dir, f"{srr_id}.fastq.gz")

    if os.path.exists(file_r2): #lecturas pair end
        return (file_r1, file_r2)
    elif os.path.exists(file_r1):
        return (file_r1, None)
    elif os.path.exists(file_s):
        return (file_s, None)
    else:
        return (None, None)

#Concatenar archivos FASTQ y comprimir 
def concatenate(output_dir, files_r1, files_r2=None, sample_name="sample"):

    for f in files_r1:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Archivo faltante: {f}")

    out_r1 = os.path.join(output_dir, f"{sample_name}_1.fastq.gz")
    if files_r2:
        out_r2 = os.path.join(output_dir, f"{sample_name}_2.fastq.gz")
    else: 
        None

    with open(out_r1, "wb") as fout:
        p1 = subprocess.Popen(["zcat"] + files_r1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
        p1.stdout.close()
        p2.communicate()


    if files_r2:
        for f in files_r2:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Archivo faltante: {f}")
            
        with open(out_r2, "wb") as fout:
            p1 = subprocess.Popen(["zcat"] + files_r2, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
            p1.stdout.close()
            p2.communicate()

#Comprueba que el archivo FASTQ esté comprimido y extension correcta
def verify_fastq(file_path):
    if not file_path.endswith((".fastq.gz", ".fastq")):
        print(f" {file_path} no es un archivo FASTQ valido.")
        return False
    return True
