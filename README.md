# Nombre del proyecto 
## High-throughput transcriptome sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth


Fecha: 29/09/2025

## Descripción del proyecto 

El presente proyecto tiene como propósito reproducir y analizar, mediante herramientas bioinformáticas, los resultados reportados en el estudio:

>Vichi J., Salazar E., Jiménez Jacinto V., Olvera Rodríguez L., Grande R., Dantán-González E., Morett E., Hernández-Mendoza A. (2021). High-throughput transcriptome >sequencing and comparative analysis of Escherichia coli and Schizosaccharomyces pombe in respiratory and fermentative growth. PLOS ONE. >https://doi.org/10.1371/journal.pone.0248513

En dicho trabajo, se llevó a cabo la secuenciación y análisis comparativo del transcriptoma de Escherichia coli (procariota) y Schizosaccharomyces pombe (eucariota unicelular) bajo condiciones de crecimiento respiratorias y fermentativas. A partir de estos datos, se identificaron genes diferencialmente expresados y se exploraron las similitudes y diferencias en los procesos metabólicos y regulatorios entre ambos organismos.

En este proyecto se emplearán directamente los datos proporcionados por los autores para llevar a cabo el pipeline de análisis bioinformático.

## Objetivo

Analizar los datos transcriptómicos de E. coli y S. pombe bajo condiciones de respiración y fermentación con el fin de identificar genes diferencialmente expresados, explorar su anotación funcional y realizar una comparación entre organismos que permita resaltar procesos conservados y divergentes.

## Planteamiento del problema

Se busca analizar los genes diferencialmente expresados en Escherichia Coli y de Schizosaccharomyces pombe mediante secuenciación transcriptómica

## Calendario de trabajo

[Definir de manera general la actividades que se requerirán para el proyecto. Por ejemplo:]

| Actividad | Fecha   | Responsable  | Entregable |
|----------|----------|----------|----------|
| Descripción de proyecto    | Reunirnos Martes y viernes de 10 a 12  | Miryam Zamora, Alondra Márquez  | Documento markdown-README.md |
| Especificación de requisitos    | septiembre   | Alo y Miryam   | Documento markdown-README.md   |
| Análisis y diseño   | 20 septiembre  |
| Construcción   | octubre, noviembre  |  Miryam y Alo    | Scripts |
| Pruebas   | noviembre  |  Miryam y Alo   | Documento markdown |
| Reporte de resultados  | noviembre  |  Alo y Miryam   | Documentos markdown |
| Presentación del proyecto   | diciembre  |  Alo y Miryam   | repositorio GitHub (release)|



## Metodología
[Descripción general de los pasos a realizar para el proyecto, por ejemplo:]

**Preguntas de investigación**



Paso 1: Localización de fuente de datos  
Paso 2: Descarga de archivos de datos  
Paso 3: Inspección de datos  
Paso 4: Limpieza de datos
Paso 5: Descripción de los datos  
Paso 5: Análsis de datos depurados  
Paso 6: Obtención de resultados  



## Resultados esperados









## Especificación de Requisitos

Requisitos funcionales

[Un requisito funcional define una función del sotware. Determina que debe hacer el software: Por ejemplo: Leer números de un archivo dado, Calcular la suma de todos los números leídos del archivo, Producir un mensaje de error si el archivo no existe, etc.]


Requisitos no funcionales

[Estos requisitos se ocupan de aspectos como la seguridad, el rendimiento, la facilidad de uso, la fiabilidad y la escalabilidad. Por ejemplo: El script deberá estar escrito en Python, El tiempo de respuesta debe ser rápido, incluso con archivos de gran tamaño, La entrada del archivo debe ser flexible (i.e. se acepta a través de la línea de comandos), etc.]




## Análisis y Diseño



Para resolver este problema, se utilizarán varias funciones incorporadas en Python, así como el manejo de excepciones para la validación de datos y archivo. A continuación, se muestra un pseudocódigo simple para ilustrar la lógica básica del script:

```
Función principal (Suma_Numero):
    Intentar:
        datos_archivo = Obtener_Datos_Archivo(ruta_archivo)
        numeros = Validar_Datos(datos_archivo)
        resultado = Calcular_Suma(numeros)
        Imprimir_Resultado(numeros, resultado)
    Atrapar cualquier excepción como error:
        Imprimir el error

Función Obtener_Datos_Archivo(ruta_archivo):
    Si la ruta del archivo no existe:
        Levantar un error de "archivo no encontrado"
    Leer y retornar las líneas del archivo

Función Validar_Datos(data):
    Intentar:
        Convertir todos los datos a formato flotante y retornar como una lista
    Atrapar ValorError:
        Levanta un error de "¡Introducir solamente números de base 10!"

Función Calcular_Suma(numeros):
    Retornar la suma de los números

Función Imprimir_Resultado(numeros, resultado):
    Imprimir los números y su suma total en el formato especificado
```

El formato de los datos de entrada será simplemente un archivo, con un número por línea. Los números pueden estar en formato entero o decimal. La salida será una línea de texto que muestra los números sumados y la suma total, en el formato: n1 + n2 + n3 + ... = suma. Los mensajes de error se imprimirán en la consola.


#### Caso de uso: Sumar Números

```
         +---------------+
         |   Usuario     |
         +-------+-------+
                 |
                 | 1. Proporciona archivo de entrada
                 v
         +-------+-------+
         |   Sumador de  |
         |   Números en  |
         |   Archivo     |
         | (Sistema)     |
         +---------------+
```

- **Actor**: Usuario
- **Descripción**: El actor proporciona un archivo de entrada con números a sumar. El sistema valida el archivo y los datos de entrada, calcula la suma de los números y muestra el resultado.
- **Flujo principal**:

	1. El actor inicia el sistema proporcionando el archivo de entrada con los números a sumar.
	2. El sistema valida el archivo y los datos de entrada.
	3. El sistema calcula la suma de los números.
	4. El sistema muestra el resultado.
	
- **Flujos alternativos**:
	- Si el archivo proporcionado no existe
		1. El sistema muestra un mensaje de error diciendo que el archivo no se encuentra.
	- Si los datos de entrada no son números en base 10
		1. El sistema muestra un mensaje de error diciendo que se deben introducir números en base 10.
