#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess

pathAnalisi = "/home/ffuster/panalisi/resultats"

def getFASTQnames(path) :
    files2copy = []
    print("INFO: Recolling el nom dels FASTQ des de {}".format(path))
    for root, dirs, files in os.walk(path) :
        for fic in files :
            aux, extension = os.path.splitext(fic)
            name, extension2 = os.path.splitext(aux)
            #Collect the FASTQ files that will be copied to "tanda" folder
            if extension2 == ".fastq" :
                pt = "{}/{}".format(root, fic)
                files2copy.append(pt)
    print("INFO: {} arxius trobats".format(len(files2copy)))

    return files2copy

def getTanda() :
    nums = []
    for root, dirs, files in os.walk(pathAnalisi) :
        for d in dirs :
            if d.startswith("tanda") :
                aux = int(d[5:])
                nums.append(aux)
        break
    sig = max(nums) + 1
    return sig

def writeBash(tanda, fastqs) :
    carpetas = []
    os.chdir(pathAnalisi)
    arx = "logTanda{}.sh".format(tanda)
    print("INFO: Creant el bash per la tanda {}".format(arx))
    with open(arx, "w") as fi :
        fi.write("#!/bin/bash\n\n")
        #Escriure els paths dels arxius usats com a referencia pels programes inclosos en la pipeline
        fi.write("#Referencies usades en aquest analisi\n")
        fi.write("ref=$HOME/panalisi/referencies/gatkHg19.fa\n")
        fi.write("mani=$HOME/panalisi/resultats/tanda28/manifest.bed\n")
        fi.write("indels=$HOME/panalisi/referencies/gold_indels.vcf\n")
        fi.write("sites=$HOME/panalisi/referencies/dbsnp_138.hg19.vcf\n")
        fi.write("gens=$HOME/panalisi/resultats/gensAestudi.txt\n")
        fi.write("manifest=$HOME/panalisi/resultats/manifest.bed\n\n")
        #Crear la funcio que copia les dades (FASTQs) en la carpeta de l'analisi
        fi.write("function copiar {\n")
        fi.write("\tcd $HOME/panalisi/resultats\n")
        fi.write("\tmkdir tanda{tanda} ; cd tanda{tanda}\n\n".format(tanda = tanda))
        fi.write("\techo -e \"################################\\n\\tCopiant dades\\n################################\\n\"\n")
        for f in fastqs :
            patx = f.replace(" ", "\\ ")
            fi.write("\trsync -aP {} .\n".format(patx))

        fi.write("\tmv ../logTanda{}.sh .\n".format(tanda))
        fi.write("}\n\n")
        ##TODO Una vegada estiguen els sripts creats per totes les etapes de la pipeline, modificar els parametres corresponents
        fi.write("function analisi {\n")
        fi.write("\tforward=$1\n\treverse=$2\n\treadgroup=$3\n\talias=$4\n")
        fi.write("\tmkdir $alias\n")
        fi.write("\tcd $alias\n")
        fi.write("\t$HOME/anpanmds/fastqc.sh ../$forward ../$reverse\n")
        fi.write("\t$HOME/anpanmds/align.sh $reference ../$forward ../$reverse $readgroup\n")
        fi.write("\t$HOME/anpanmds/postAlign.sh recalibrate bwaAlign/bwa.sort.bam\n")
        fi.write("coverage")
        fi.write("\t$HOME/anpanmds/variantCalling.sh $manifest bwaAlign/bwa.recalibrate.bam\n")
        fi.write("\T$HOME/anpanmds/variantAnnotation.sh variants.vcf")
        fi.write("}\n\n")

        ## TODO: Canviar el format de la pipeline
        # Crear una estructura similar a l'exemple de l'editor de textos
        # Afegir a tots els scripts la capçalera "script utilitzat a partir de la tanda33. Els anteriors estan en la carpeta panalis/scripts"


"""
        fi.write("echo -e \"################################\\n\\tComençant analisi\\n################################\\n\"\n")
        for f in fastqs :
            i = f.split('/')[-1]
            try :
                id = int(i.split('_')[0])
                samp = i.split('_')[1]
                rg = "\"@RG\\tID:{}-IDH\\tSM:{}\\tPL:ILLUMINA\"".format(id, samp)
                outputF = "{}-IDH".format(id)
                f1 = "{}_{}_L001_R1_001.fastq.gz".format(id, samp)
                f2 = "{}_{}_L001_R2_001.fastq.gz".format(id, samp)
            except ValueError, e :
                id = i.split('_')[0]
                samp = i.split('_')[1]
                rg = "\"@RG\\tID:{}\\tSM:{}\\tPL:ILLUMINA\"".format(id, samp)
                outputF = id
                f1 = "{}_{}_L001_R1_001.fastq.gz".format(id, samp)
                f2 = "{}_{}_L001_R2_001.fastq.gz".format(id, samp)

            if not outputF in carpetas :
                fi.write("./analizar.sh {} {} {} {}\n".format(f1, f2, rg, outputF))
                carpetas.append(outputF)

        fi.write("\necho -e \"################################\\n\\tCopiant resultats en servidor\\n################################\\n\"\n")
        for c in carpetas :
            fi.write("rsync -aP {0}/bwaAlign/bwa.realigned.ba* {0}/variantCalling/{0}.xlsx /home/labs/solelab/ffuster2/Desktop/sequence/PANELLS/FSole/Myeloid\ Panel_Agilent_2017/Analisi_Fco/{0}\n".format(c))

    subprocess.call("chmod u+x logTanda{}.sh".format(tanda), shell=True)
    print "INFO: Bash created, named as {}".format(arx)
    print "INFO: Executing {} as a subprocess".format(arx)
    com = "./{}".format(arx)
    pr = subprocess.Popen(com, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logIn, logErr = pr.communicate()
    if pr.returncode == 0 :
        print "INFO: Analysis finished. No errors found"
    else :
        print "ERROR: Analysis ended with errors. Check logs"

    with open("tanda{}.log".format(tanda), "w") as fi :
        fi.write(logIn)
    with open("tanda{}.err".format(tanda), "w") as fi :
        fi.write(logErr)

    print "INFO: Logs written as tanda{0}.log & tanda{0}.err".format(tanda)
"""


def main() :
    path = input("\n\tREQUEST: Dona'm el path absolut on estan els FASTQ de la tanda a analitzar (Consell: Fes un `cd` a la carpeta i executar `pwd`)\n")
    if os.path.isdir(path) :
        files2copy = getFASTQnames(path)
        fold = getTanda()
        print("INFO: L'analisi es guardara en {ruta}/tanda{num}".format(ruta = pathAnalisi, num = fold))
        writeBash(fold, files2copy)
    else :
        print("ERROR: No he trobat el directori {}".format(path))

if __name__ == "__main__" :
    main()

"""

def writeBash(tanda, fastqs) :
    carpetas = []
    os.chdir(anPath)
    arx = "logTanda{}.sh".format(tanda)
    print "INFO: Creating bash {}".format(arx)
    with open(arx, "w") as fi :
        fi.write("#!/bin/bash\n\n")
        fi.write("cd $HOME/panalisi/resultats\n")
        fi.write("mkdir tanda{0} ; cd tanda{0}\n\n".format(tanda))
        fi.write("echo -e \"################################\\n\\tCopiant dades\\n################################\\n\"\n")
        for f in fastqs :
            patx = f.replace(" ", "\\ ")
            fi.write("rsync -aP {} .\n".format(patx))

        fi.write("\necho -e \"################################\\n\\tCopiant arxius necessaris per analisi\\n################################\\n\"\n")
        fi.write("rsync -aP ../tanda{0}/analizar.sh ../tanda{0}/gensAestudi.txt ../tanda{0}/manifest.bed .\n".format(tanda-1))
        fi.write("mv ../logTanda{}.sh .\n\n".format(tanda))

        fi.write("echo -e \"################################\\n\\tComençant analisi\\n################################\\n\"\n")
        for f in fastqs :
            i = f.split('/')[-1]
            try :
                id = int(i.split('_')[0])
                samp = i.split('_')[1]
                rg = "\"@RG\\tID:{}-IDH\\tSM:{}\\tPL:ILLUMINA\"".format(id, samp)
                outputF = "{}-IDH".format(id)
                f1 = "{}_{}_L001_R1_001.fastq.gz".format(id, samp)
                f2 = "{}_{}_L001_R2_001.fastq.gz".format(id, samp)
            except ValueError, e :
                id = i.split('_')[0]
                samp = i.split('_')[1]
                rg = "\"@RG\\tID:{}\\tSM:{}\\tPL:ILLUMINA\"".format(id, samp)
                outputF = id
                f1 = "{}_{}_L001_R1_001.fastq.gz".format(id, samp)
                f2 = "{}_{}_L001_R2_001.fastq.gz".format(id, samp)

            if not outputF in carpetas :
                fi.write("./analizar.sh {} {} {} {}\n".format(f1, f2, rg, outputF))
                carpetas.append(outputF)

        fi.write("\necho -e \"################################\\n\\tCopiant resultats en servidor\\n################################\\n\"\n")
        for c in carpetas :
            fi.write("rsync -aP {0}/bwaAlign/bwa.realigned.ba* {0}/variantCalling/{0}.xlsx /home/labs/solelab/ffuster2/Desktop/sequence/PANELLS/FSole/Myeloid\ Panel_Agilent_2017/Analisi_Fco/{0}\n".format(c))

    subprocess.call("chmod u+x logTanda{}.sh".format(tanda), shell=True)
    print "INFO: Bash created, named as {}".format(arx)
    print "INFO: Executing {} as a subprocess".format(arx)
    com = "./{}".format(arx)
    pr = subprocess.Popen(com, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logIn, logErr = pr.communicate()
    if pr.returncode == 0 :
        print "INFO: Analysis finished. No errors found"
    else :
        print "ERROR: Analysis ended with errors. Check logs"

    with open("tanda{}.log".format(tanda), "w") as fi :
        fi.write(logIn)
    with open("tanda{}.err".format(tanda), "w") as fi :
        fi.write(logErr)

    print "INFO: Logs written as tanda{0}.log & tanda{0}.err".format(tanda)

def main() :
    #Ask for the folder
    inp = raw_input("Check in {} the new folder to analyse.\nWrite the folder (relative path) where the FASTQs are: ".format(absPath))
    fullP = "{}/{}".format(absPath, inp)

    #TODO fer flags per cadascun dels programes
"""
