# Analisis Utom&aacute;tico de Paneles (sic)
## Analizar paneles SMD autom&aacute;ticamente

Procesa las muestras de paneles secuenciadas en genomica autom&aacute;ticamente

## Uso

```bash
python3 main.py
```
El programa pedir&aacute; una ruta absoluta donde est&aacute;n los FASTQ reportados por el MiSeq. Con introducir ese dato basta para hacer todo el an&aacute;lisis

## &iquest;C&oacute;mo hago para modificar los programas que se ejecutan en los an&aacute;lisis?

Para modificar los programas, o las versiones de programas que se est&aacute;n ejecutando en el panel, hay que modificar las constantes que aparecen al principio de [creaLog.py]{../master/creaLog.py}.

## Otras funcionalidades
### Limpiar una tanda para archivarla

```bash
python3 cleaner.py
```
