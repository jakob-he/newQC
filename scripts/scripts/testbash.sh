#!/bin/bash

for filename in /home/jakobh/Dokumente/PEDIA/PEDIA-workflow-new_download/1_qualityCheck/json/mapped/*.json; do
    python3 newQC.py --jsonfile "$filename" --jsoncurratedfolder "./../json/currated" --errordict "./../hgvs_errordict.json" --logfile "./../QC_progress.xls" --debugfolder "./../json/debug" --vcffolder "./../json/vcf"
done
