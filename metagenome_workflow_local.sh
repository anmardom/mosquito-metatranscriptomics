#!/bin/bash

# Workflow para analizar cualquier fichero transcriptomica hasta R.

# Ayuda para cuando no se usen variables correctamente
usage(){
echo -e "\nPlease, specify SRA code (Ej: SRR12064742), Kraken2 database path and hisat2 index path to execute properly. Specify if it's single end with -se (default for paired end)"
echo -e "\nmetagen_hpc.sh [-sra code] [-dbk2 /path/db] [-index /path/hisat2-index] {-se}"
echo -e "\nRecuerda ajustar paths de variables a los programas. Probar SRA del SRR12064736 al SRR12064743 (son de Pepe)"
exit 1
}


# Paso 0: Prepara variables y directorios

se="FALSE"

# Parsea variables
while [[ $# -gt 0 ]]; do
    case $1 in 
        -sra) sra="$2"; shift;;					# Nombre fichero SRA
        -se) se="TRUE";; 					# Si se indica,se tratara como SE
        -dbk2) db=$2; shift;;					# Database Kraken2
        -index) index=$2; shift;;				# Indice A. Gambiae alineamiento Hisat2
	-h|--help) usage;;
        *) echo "Opcion invalida: $1"; usage;;						
    esac
    shift
done


# Programas
fqd="programas/sratoolkit.3.0.10-ubuntu64/bin/fasterq-dump"	# Fasterq-dump
fastQC="/home/alumno/BioInfo-Programas/FastQC/fastqc"
trim="/home/alumno/BioInfo-Programas/Trimmomatic-0.39/trimmomatic-0.39.jar"
hs2="programas/hisat2-2.2.1/hisat2"				# Hisat2
k2="programas/kraken2-2.1.3/kraken2"				# Kraken2
b2="programas/Bracken-2.7/bracken"				# Bracken

res="results/${sra}"						# Directorio de resultados
mkdir "${res}"
index="reference_genomes_db/index_mosquito/agambiae_index"	# Genoma A. Gambiae referencia alineamiento
db="databases/minikraken_8GB_20200312/"				# Database referencia (local, mala)
adapters="/home/alumno/BioInfo-Programas/Trimmomatic-0.39/adapters/TruSeq2-PE.fa"



# Paso 1: Obtener fichero SRR (utilizando SRAtoolkit - fasterq-dump junto al nombre del fichero)

echo -e "\nFasterq-dumping: ${sra}" 
$fqd $sra -O $res  

echo "Fichero ${sra} descargado" 



# Paso 2: Analisis de Calidad FastQC (Mostrar resultados y preguntar si continuar o trimmear)

if [ "$se" = "TRUE" ]; then
    srr="${res}/${sra}.fastq"					 # Fichero SE
    
    ${fastQC} ${srr} 						 
    xdg-open "${res}/${sra}_fastqc.html" 			 # Abre html en firefox
    #read -p "Continue? (n: )" choice
    #if [ "$choice" = "n" ]; then
    #    exit
    #fi

    # Elimina Adaptadores
    echo -e "\nTrimming SE: ${sra}" 
    java -jar "${trim}" SE -phred33 -threads 4 "${srr}" "${res}/SE.fastq" ILLUMINACLIP:/home/alumno/BioInfo-Programas/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 

    echo "Trim completado" 

    # Actualiza nombre de ficheros
    rm $srr
    srr="${res}/SE.fastq"
    ${fastQC} ${srr} 

else
    srr1="${res}/${sra}_1.fastq"				 # Fichero PE 1
    srr2="${res}/${sra}_2.fastq"				 # Fichero PE 2
    
    ${fastQC} ${srr1} 
    ${fastQC} ${srr2} 
    #xdg-open "${res}/${sra}_1_fastqc.html" 		 	 # Abre html 1 en firefox
    #xdg-open "${res}/${sra}_2_fastqc.html" 		 	 # Abre html 2 en firefox
    #read -p "Continue? (n: )" choice
    #if [ "$choice" = "n" ]; then
    #    exit
    #fi
    
    # Elimina adaptadores
    echo -e "\nTrimming PE: ${sra}" 
    java -jar "${trim}" PE -phred33 -threads 4 "${srr1}" "${srr2}" "${res}/1_PE.fastq" "${res}/1_uPE.fastq" "${res}/2_PE.fastq" "${res}/2_uPE.fastq" ILLUMINACLIP:${adapters}:2:30:10 

    echo "Trim completado" 

    # Actualiza nombre de ficheros
    rm $srr1
    rm $srr2
    srr1="${res}/1_PE.fastq"
    srr2="${res}/2_PE.fastq"
    ${fastQC} ${srr1} 
    ${fastQC} ${srr2} 
    
fi



# Paso 3: Alineamiento Hisat2 con genomas de A. Gambiae

echo -e "\nAlineando con el genoma de referencia..." 

    # Si es SE:

if [ "$se" = "TRUE" ]; then
    ${hs2} -x ${index} -p 4 -U ${srr} | samtools view -b > "${res}/alineamiento_${sra}.bam" 
    rm $srr

    # Si es PE: 
else
    ${hs2} -x ${index} -p 4 -1 ${srr1} -2 ${srr2} | samtools view -b > "${res}/alineamiento_${sra}.bam" 
    rm $srr1
    rm $srr2
fi

    # Alternativa (descartar, saltaría paso Quality Control):
#${hs2} -x ${index} -p 4 -1 --sra-acc "$srr" | samtools view -b > "${res}/alineamiento_${srr}_hisat2.bam"

echo "Alineamiento completado" 



# Paso 4: Analisis SAMtools y extraer secuencias no alineadas en formato fastq (SE)

echo -e "\nObteniendo secuencias ALINEADAS en formato bam..." 
samtools view -b -F 4 "${res}/alineamiento_${sra}.bam" -o "${res}/alligned_${sra}.bam" 

echo -e "\nObteniendo secuencias NO ALINEADAS en formato fastq.gz..." 
samtools view -b -f 4 "${res}/alineamiento_${sra}.bam" > "${res}/unalligned_${sra}.bam" 
samtools fastq "${res}/unalligned_${sra}.bam" | gzip -k > "${res}/unalligned_${sra}.fastq.gz" 
rm "${res}/alineamiento_${sra}.bam"
rm "${res}/unalligned_${sra}.bam"

echo "Secuencias extraidas" 

# Alternativa para secuencias PE
#gzip -dk unalligned_SRR12064741.fastq.gz

#samtools fastq -f 64 unalligned_SRR12064741.fastq > unalligned_SRR12064741_R1.fastq
#gzip unalligned_SRR12064741_R1.fastq

#samtools fastq -f 128 unalligned_SRR12064741.fastq > unalligned_SRR12064741_R2.fastq
#gzip unalligned_SRR12064741_R2.fastq

# Convierte en krona y compara
#programas/Krona-2.8/KronaTools/scripts/ImportTaxonomy.pl -t 2 -m 6 -o results/SRR12064740/krona_SEvsPE.html -tax programas/Krona-2.8/taxonomy/ results/SRR12064740/meta2_SRR12064740_bracken77GB_S.bracken,PE c3upo_trials/results/bracken77GB_SRR12064740_S.bracken,SE



# Paso 5: Analisis Kraken2

#echo -e "\nAnalisis taxonomico con Kraken2 y minikraken (8 GB)"
#$k2 --use-names --threads 4 --db "${db}" --output "${res}/results_unalligned_${sra}.kraken2" --report "${res}/report_unalligned_${sra}.kreport2" --gzip-compressed "${res}/unalligned_${sra}.fastq.gz"



# Paso 6: Analisis Bracken (probar filtro -t 0/5/10)

#$b2 -d "${db}" -i "${res}/report_unalligned_${sra}.kreport2" -o "${res}/bracken_unalligned_${sra}_S.bracken" -w "${res}/bracken_askraken2_unalligned_${sra}_S.kreport2" -r 150 -l S -t 10



# Paso 7: Procesar ficheros en resultado conjunto (t.sh para 2 muestras, jupyter-lab para n muestras -> IGNORAR)

#awk -F'\t' 'NR==FNR {if (FNR>1) a[$1]=$6; next} {if (FNR>1) b[$1]=$6} END { for (otu in a) {print otu, (otu in b) ? a[otu] : 0, (otu in b) ? b[otu] : 0} for (otu in b) if (!(otu in a)) print otu, 0, b[otu] }' ${res}/bracken_unalligned_${srr}_hisat2_species.bracken ${res}/bracken_unalligned_${srr}_hisat2_species.bracken > trial.bracken















