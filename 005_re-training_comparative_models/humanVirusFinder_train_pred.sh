#!/bin/bash

### Before running, please setup your environment following the git-hub (https://github.com/Fzhang1992/human-infecting_virus_finder)

mkdir output_examples/humanVirusFinder

fold="fold_0"

## preparing sequence 
# Coronaviridae.curated.csv: the fold_0 information is present in column 3 and the infectivity label in column 4.
for dataset in "train" "dev" "pred"
do

cat ../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt | \
	awk -v dataset=$dataset 'BEGIN{FS=","} NR>1{if($3 == dataset && $4 == 0) print $2}' \
	> output_examples/humanVirusFinder/${dataset}_animal.txt

cat ../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt | \
	awk -v dataset=$dataset 'BEGIN{FS=","} NR>1{if($3 == dataset && $4 == 1) print $2}' \
	> output_examples/humanVirusFinder/${dataset}_human.txt

cat ../002_Data_splitting/input_examples/Coronaviridae.curated.fasta | \
	seqkit grep -f output_examples/humanVirusFinder/${dataset}_animal.txt \
	> output_examples/humanVirusFinder/${dataset}_animal.fasta

cat ../002_Data_splitting/input_examples/Coronaviridae.curated.fasta | \
	seqkit grep -f output_examples/humanVirusFinder/${dataset}_human.txt \
	> output_examples/humanVirusFinder/${dataset}_human.fasta

done

## training
python human-infecting_virus_finder/train_mod.py \
	--infect output_examples/humanVirusFinder/train_human.fasta \
	--other output_examples/humanVirusFinder/train_animal.fasta \
	--file output_examples/humanVirusFinder/model_${fold} \
	--kmer 4


## Prediction
cat output_examples/humanVirusFinder/pred_*.fasta > output_examples/humanVirusFinder/pred.fasta

python human-infecting_virus_finder/predResult.py \
	--query output_examples/humanVirusFinder/pred_${fold}.fasta \
	--output output_examples/humanVirusFinder/pred_${fold}.txt \
	--file output_examples/humanVirusFinder/model_${fold} \
	--kmer 4 \
	--kjob 1 \
	--bjob 1