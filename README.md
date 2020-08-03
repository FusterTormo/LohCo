# Analisis aUtom&aacute;tico de Paneles
## Analizar paneles SMD autom&aacute;ticamente

Procesa las muestras de paneles secuenciadas en genomica autom&aacute;ticamente

## Uso

```bash
python3 main.py
```
El programa pedir&aacute; una ruta absoluta donde est&aacute;n los FASTQ reportados por el MiSeq. Con introducir ese dato basta para hacer todo el an&aacute;lisis, aunque se pueden a&ntilde;adir los pasos para realizar un an&aacute;lisis m&aacute;s customizado.

# FAQs

## &iquest;C&oacute;mo hago para modificar los programas que se ejecutan en los an&aacute;lisis?

La liber&iacute;a [getCommands.py](../master/getCommands.py) contiene las &oacute;rdenes para ejecutar cada uno de los programas necesarios para el an&aacute;lisis. En las constantes de esta librer&iacute;a se encuentran las versiones de todos los programas que se usan para el an&aacute;lisis.


# Otras funcionalidades
## Limpiar una tanda para archivarla

Eliminar todos los archivos no importantes generados durante el an&aacute;lisis de una tanda de paneles

### Uso

```bash
python3 cleaner.py [numero_de_tanda]
```

El programa eliminar&aacute; los archivo temporales de la tanda pasada por par&aacute;metro. Los archivos temporales son los FASTQ, los bams intermedios y los archivos anotados intermedios.

```bash
python3 cleaner.py
```

El programa pedir&aacute; el n&uacute;mero de tanda que se quiere limpiar. Una vez hecho esto, ir&aacute; a la carpeta correspondiente y eliminar&aacute; los archivos FASTQ, los bams intermedios y los vcf y anotaciones intermedias.

## Analizar una muestra por separado
## Analisis personalizado
## Reanalizar una muestra sin eliminar los datos hechos previamente
## Crear un excel a partir de un vcf de variantes

## Control de calidad de un alineamiento en particular

El control de calidad se puede ejecutar de dos formas distintas: invocando directamente al *script* de control de calidad (bamQC.py) o usando la interfaz (main.py)

### Usando bamQC

El control de calidad se debe hacer dentro de la carpeta de la muestra. Por esto, para ejecutar este *script* se debe poner en la carpeta de la muestra. Una vez all&iacute; ejecutamos el *script*

```bash
cd $MUESTRA_A_ANALIZAR
python3 $PATH_AUP/bamQC.py
```

El programa crear&aacute; un archivo de texto, en formato JSON, que contiene el n&uacute;mero de *reads* que tiene el FASTQ (suma de los *reads* reportados por FASTQC), los *reads* que tiene el bam (sumatorio de las filas que tiene el bed creado), los *reads* que est&aacute;n que caen dentro del manifest (ON target) y los que no caen dentro del manifest (OFF target).

### Usando main

Ejecutar el programa principal

```bash
python3 main.py
```

La interfaz pedir&aacute; indicar la ruta donde est&aacute; el bam en el que se quiere hacer el control de calidad. Con solo este dato, la interfaz crear&aacute; los archivos necesarios para poder hacer dicho control de calidad y reportar&aacute; la informaci&oacute;n del archivo en formato JSON en un formato m&aacute;s legible. Una vez se hayan mostrado los resultados por pantalla, la interfaz borrar&aacute; todos los archivos que ha creado.
