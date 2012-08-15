#!/bin/sh

if [ $# != 4 ]
then
    echo "Usage: amit_var_filter.sh <GATK BIN> <REF> <INPUT_DIR> <FILENAME> ";
else
	gatk=$1
	ref=$2
	out=$3
	filename=$4

	/usr/java/latest/bin/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
	-R $ref -et NO_ET -K $gatk/Hossain.Asif_mayo.edu.key -l INFO \
	-T VariantFiltration -V $out/$filename -o $out/$filename.tmp \
	--filterExpression "QD < 2.0" --filterName QDFilter \
	--filterExpression "HRun > 5" --filterName HRunFilter \
	--filterExpression "FS > 20" --filterName FSFilter

	mv $out/$filename.tmp $out/$filename
fi
