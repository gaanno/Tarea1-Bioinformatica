﻿# Tarea1-Bioinformatica

Asunto: Tarea 1 Bioinformática

Programa de alineamiento de secuencias de pares utilizando BLOSUM62 y PAM250.

- Busca las secuencias en GenBank.
- Realiza alineamiento global y local de cada par.
- Muestra score y porcentaje de identidad. 

## Tabla de contenidos
1. [Requerimientos](#requerimientos)
2. [Uso](#uso)
3. [Información API GenBank](#información-api-genbank)
4. [Bases de datos disponibles](#bases-de-datos-disponibles)


## Requerimientos

Para instalar las bibliotecas requeridas

~~~
pip3 install -r requirements.txt
~~~

## Uso
Colocar el nombre de las secuencias en pares que se desean buscar en la base de datos de GenBank en un archivo separados por espacios con formato

~~~
acc acc
acc acc
acc acc
~~~

En donde:
    acc: Corresponde al accession number

Ejecutar programa con 

~~~python
python3 main.py -file nucleotide.txt -type nucleotide -format fasta


python3 main.py -file protein.txt -type protein -format fasta
~~~

En donde:
1. -file: (Obligatorio) Archivo contenedor de las acc.
2. -type: (Obligatorio) Base de datos a consultar.
3. -format: (Opcional)  Formato que se desea usar como fasta, fa, gb, etc.

## Información API GenBank

Link de la API de GenBank https://www.ncbi.nlm.nih.gov/home/develop/api/

## Bases de datos disponibles

Para consultar las bases de datos disponibles https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi
