#!/bin/bash

### Before running, please setup your environment following the git-hub (https://github.com/flomock/vidhop)

mkdir output_examples/vidhop

fold="fold_0"

## preparing sequence 
# Coronaviridae.curated.csv: the fold_0 information is present in column 3 and the infectivity label in column 4.
for dataset in "train" "dev" "pred"
do

cat ../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt | \
	awk -v dataset=$dataset 'BEGIN{FS=","} NR>1{if($3 == dataset && $4 == 0) print $2}' \
	> output_examples/vidhop/${dataset}_animal.txt

cat ../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt | \
	awk -v dataset=$dataset 'BEGIN{FS=","} NR>1{if($3 == dataset && $4 == 1) print $2}' \
	> output_examples/vidhop/${dataset}_human.txt

cat ../002_Data_splitting/input_examples/Coronaviridae.curated.fasta | \
	seqkit grep -f output_examples/vidhop/${dataset}_animal.txt  \
	> output_examples/vidhop/${dataset}_animal.fasta

cat ../002_Data_splitting/input_examples/Coronaviridae.curated.fasta | \
	seqkit grep -f output_examples/vidhop/${dataset}_human.txt  \
	> output_examples/vidhop/${dataset}_human.fasta

cat output_examples/vidhop/${dataset}_animal.fasta \
	output_examples/vidhop/${dataset}_human.fasta \
	> output_examples/vidhop/${dataset}.fasta

# convert to vidhop inputs
python \
	VIDHOP/vidhop_input.py \
	output_examples/vidhop/${dataset}.fasta \
	../003_Data_predictability/input_examples/known_virus_final_dataset_Coronaviridae.txt \
	output_examples/vidhop/${dataset}_sequences.txt \
	output_examples/vidhop/${dataset}_hosts.txt

done

## make dataset (train)
vidhop make_dataset \
	-x output_examples/vidhop/train_sequences.txt \
	-y output_examples/vidhop/train_hosts.txt \
	--repeated_undersampling \
	--val_split_size 0 \
	--test_split_size 0 \
	-o output_examples/vidhop/train_data/

## make dataset (eval)
vidhop \
	make_dataset \
	-x output_examples/vidhop/dev_sequences.txt \
	-y output_examples/vidhop/dev_hosts.txt \
	--repeated_undersampling \
	--val_split_size 1 \
	--test_split_size 0 \
	-o output_examples/vidhop/train_data/

## make dataset (test)
vidhop \
	make_dataset \
	-x output_examples/vidhop/pred_sequences.txt \
	-y output_examples/vidhop/pred_hosts.txt \
	--repeated_undersampling \
	--val_split_size 0 \
	--test_split_size 1 \
	-o output_examples/vidhop/train_data/

## training
vidhop training \
	-i output_examples/vidhop/train_data \
	-o output_examples/vidhop/model_${fold} \
	--name "Coronaviridae" \
	--architecture 0 \
	--epochs 1500  --early_stopping


## prediction
vidhop predict \
	--input output_examples/vidhop/pred.fasta \
	--virus output_examples/vidhop/model_${fold}/model_best_auc_Coronaviridae.model \
	--outpath output_examples/vidhop/pred_output_${fold}