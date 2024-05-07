#!/bin/bash

### Before running, please setup your environment following the git-hub (https://gitlab.com/dacs-hpi/deepac-vir)

mkdir output_examples/deepac

fold="fold_0"

## preparing sequence 
# Coronaviridae.curated.csv: the fold_0 information is present in column 3 and the infectivity label in column 4.
for dataset in "train" "dev" "pred"
do

cat ../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt | \
	awk -v dataset=$dataset 'BEGIN{FS=","} NR>1{if($3 == dataset && $4 == 0) print $2}' \
	> output_examples/deepac/${dataset}_animal.txt

cat ../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt | \
	awk -v dataset=$dataset 'BEGIN{FS=","} NR>1{if($3 == dataset && $4 == 1) print $2}' \
	> output_examples/deepac/${dataset}_human.txt

cat ../002_Data_splitting/input_examples/Coronaviridae.curated.fasta | \
	seqkit grep -f output_examples/deepac/${dataset}_animal.txt | \
	seqkit sliding -s 125 -W 250 --greedy | \
	seqkit seq -m 100 \
	> output_examples/deepac/${dataset}_animal.fasta

cat ../002_Data_splitting/input_examples/Coronaviridae.curated.fasta | \
	seqkit grep -f output_examples/deepac/${dataset}_human.txt | \
	seqkit sliding -s 125 -W 250 --greedy | \
	seqkit seq -m 100 \
	> output_examples/deepac/${dataset}_human.fasta

deepac preproc DeePac_vir/preproc_${dataset}.ini

done

## training
class0=$(cat output_examples/deepac/train_animal.fasta | grep ">" | wc -l)
class1=$(cat output_examples/deepac/train_human.fasta | grep ">" | wc -l)

cat DeePac_vir/nn-img-rapid-cnn_todai.ini | \
	sed -e "s/{class0}/$class0/g" -e "s/{class1}/$class1/g"  \
	> DeePac_vir/nn-img-rapid-cnn_class.ini
	
deepac train -c DeePac_vir/nn-img-rapid-cnn_class.ini

## prediction
model=$(ls logs/img-rapid-cnn-logs/*.h5 | grep -oP 'e\d{3}' | sort -n | tail -1 | xargs -I {} find logs/img-rapid-cnn-logs -name "*{}*.h5")

deepac predict \
	-a output_examples/deepac/pred_data_${fold}.npy \
	-c  ${model}  \
	-o output_examples/deepac/pred_${fold}.npy