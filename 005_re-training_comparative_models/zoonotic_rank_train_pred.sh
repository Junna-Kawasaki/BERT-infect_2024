#!/bin/bash

### Before running, please setup your environment following the git-hub (https://gitlab.com/dacs-hpi/deepac-vir)
## please prepare "Sequence" folder following the git-hub.

mkdir output_examples/zoonotic_rank

train_virus="Coronaviridae"
fold="fold_0"

## extracting genomic features
Rscript \
	zoonotic_rank/Zoonotic_rank_006_CalculateGenomicFeatures.R \
	zoonotic_rank/FinalData_Cleaned_${fold}.csv \
	zoonotic_rank/Sequences \
	output_examples

## training
Rscript \
	script/Zoonotic_rank_008a_TrainAndValidate_FeatureSelection_Top_0808.R \
	output_examples

## prediction
Rscript \
	script/Zoonotic_rank_012_PredNovel.R \
	output_examples \
	zoonotic_rank/FinalData_Cleaned_${fold}.csv \
	zoonotic_rank/Sequence/AB257344.gb \ 
	test_metadata.csv