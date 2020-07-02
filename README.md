# Analisis Utom&aacute;tico (sic) de Paneles
## Analizar paneles SMD autom&aacute;ticamente

Procesa las muestras de paneles secuenciadas en genomica autom&aacute;ticamente

## Uso

```bash
python3 main.py
```
El programa pedir&aacute; una ruta absoluta donde est&aacute;n los FASTQ reportados por el MiSeq. Con introducir ese dato basta para hacer todo el an&aacute;lisis, aunque se pueden a&ntilde;adir los pasos para realizar un an&aacute;lisis m&aacute;s customizado.

# FAQs

## &iquest;C&oacute;mo hago para modificar los programas que se ejecutan en los an&aacute;lisis?

La liber&iacute;a [getCommands.py]{../master/getCommands.py} contiene las &oacute;rdenes para ejecutar cada uno de los programas necesarios para el an&aacute;lisis. En las constantes de esta librer&iacute; se encuentran las versiones de todos los programas que se usan para el an&aacute;lisis.


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
## Analisis customizado
## Reanalizar una muestra sin eliminar los datos hechos previamente
## Crear un excel a partir de un vcf de variantes
## Control de calidad de un alineamiento en particular
