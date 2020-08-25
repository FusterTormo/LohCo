# Analisis aUtom&aacute;tico de Paneles
## Analizar paneles SMD autom&aacute;ticamente

Procesa las muestras de paneles secuenciadas en genomica autom&aacute;ticamente

## Uso

```bash
$PATH_AUP/main.py
```
La interfaz pedir&aacute; una ruta absoluta donde est&aacute;n los FASTQ reportados por el secuenciador. Con introducir ese dato basta para hacer todo el an&aacute;lisis, aunque se pueden a&ntilde;adir los pasos para realizar un an&aacute;lisis m&aacute;s customizado.

## Otras funcionalidades
### Limpiar una tanda para archivarla

Eliminar todos los archivos no importantes generados durante el an&aacute;lisis de una tanda de paneles

#### Uso

```bash
$PATH_AUP/cleaner.py [numero_de_tanda]
```

El programa eliminar&aacute; los archivo temporales de la tanda pasada por par&aacute;metro. Los archivos temporales son los FASTQ, los bams intermedios y los archivos anotados intermedios.

```bash
$PATH_AUP/cleaner.py
```

El programa pedir&aacute; el n&uacute;mero de tanda que se quiere limpiar. Una vez hecho esto, ir&aacute; a la carpeta correspondiente y eliminar&aacute; los archivos FASTQ, los bams intermedios y los vcf y anotaciones intermedias.

### Anotar un manifest

En caso de tener un manifest en el que no se sabe a qu&eacute; genes pertenece cada regi&oacute;n se puede usar el *main* (la interfaz) para anotar dicho manifest. El *script*, adem&aacute;, crear&aacute; el archivo gensAestudi.txt que se usa durante la fase de c&aacute;lculo de *coverage*.

#### Uso

```bash
$PATH_AUP/main.py
```

&rarr; Opci&oacute;n n&uacute;mero 7

### Analizar una muestra por separado

Analiza una muestra que no ha sido analizada previamente.

#### Uso

```bash
$PATH_AUP/main.py
```

&rarr; Opci&oacute;n n&uacute;mero 3. La interfaz pedir&aacute; la ruta absoluta de la carpeta donde se encuentran los FASTQ que se quieren analizar. Se crear&aacute; el log del an&aacute;lisis y, posteriormente se ejecutar&aacute; dicho an&aacute;lisis.

## Analisis personalizado

Usando la interfaz (main.py), el programa preguntar&aacute; qu&eacute; pasos realizar. Solo  hay que ir contestando 's' o 'n' a cada una de las preguntas que hace el programa. Una vez decidido el an&aacute;lisis, el programa preguntar&aacute; por la ruta donde est&aacute;n las muestras a analizar. Autom&aacute;ticamente se crear&aacute; el log para analizar las muestras con los pasos que se han especificado.

## Reanalizar una muestra sin eliminar los datos hechos previamente

Usando la interfaz (main.py), el programa pedir&aacute; el identificador de la muestra que se quiere reanalizar. El programa buscar&aacute; la carpeta de la tanda en la que se analiz&oacute; la muestra. Si se encuentra dicha carpeta, el programa preguntar&aacute; por el an&aacute;lisis personalizado que se quiere hacer y crear&aacute; el log con todas las instrucciones necesarias para hacer dicho an&aacute;lisis. Como el an&aacute;lisis se har&aacute; en todas las muestras que est&eacute;n en la carpeta, si hay alguna muestra que no se quiera reanalizar, se recomienda eliminar dichas l&iacute;neas del log de an&aacute;lisis.

## Crear un excel a partir de un vcf de variantes

Usando la interfaz (main.py), pedir&aacute; la ruta absoluta del vcf que se quiere anotar. Adem&aacute;s pedir&aacute; que se le d&eacute; un nombre identificativo a la muestra que se va a anotar. Finalmente, se pedir&aacute; el genoma de referencia que se ha usado para hacer el an&aacute;lisis previo.

## Control de calidad de un alineamiento en particular

El control de calidad se puede ejecutar de dos formas distintas: invocando directamente al *script* de control de calidad (bamQC.py) o usando la interfaz (main.py)

### Usando bamQC

El control de calidad se debe hacer dentro de la carpeta de la muestra. Por esto, para ejecutar este *script* se debe poner en la carpeta de la muestra. Una vez all&iacute; ejecutamos el *script*

```bash
cd $MUESTRA_A_ANALIZAR
$PATH_AUP/bamQC.py
```

El programa crear&aacute; un archivo de texto, en formato JSON, que contiene el n&uacute;mero de *reads* que tiene el FASTQ (suma de los *reads* reportados por FASTQC), los *reads* que tiene el bam (sumatorio de las filas que tiene el bed creado), los *reads* que est&aacute;n dentro del manifest (ON target) y los que no est&aacute;n dentro del manifest (OFF target).

### Usando main

Ejecutar el programa principal

```bash
$PATH_AUP/main.py
```

La interfaz pedir&aacute; indicar la ruta donde est&aacute; el bam en el que se quiere hacer el control de calidad. La interfaz crear&aacute; los archivos necesarios para poder hacer dicho control de calidad y reportar&aacute; la informaci&oacute;n del archivo en formato JSON en un formato m&aacute;s legible. Una vez se hayan mostrado los resultados por pantalla, la interfaz borrar&aacute; todos los archivos que ha creado.

# FAQs

## &iquest;C&oacute;mo hago para modificar los programas que se ejecutan en los an&aacute;lisis?

La liber&iacute;a [getCommands.py](../master/getCommands.py) contiene las &oacute;rdenes para ejecutar cada uno de los programas necesarios para el an&aacute;lisis. En las constantes de esta librer&iacute;a se encuentran las versiones de todos los programas que se usan para el an&aacute;lisis.


## El manifest del panel ha cambiado. &iquest;C&oacute;mo cambio el manifest de la pipeline?

Existen dos opciones:

1. Crear el *log* de an&aacute;lisis (usando los comandos descritos al inicio de este README) y cambiar manualmente la ruta del manifest. Dicha ruta est&aacute; guardada como una constante en el *log* del an&aacute;lisis.
2. Cambiar la ruta de la constante manifest en [bamQC](../master/bamQC.py), y [creaLog](../master/creaLog.py)
