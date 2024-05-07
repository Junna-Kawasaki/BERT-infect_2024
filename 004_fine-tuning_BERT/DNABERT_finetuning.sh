#!/bin/bash

export tokens=250
export LR=2e-4
export KMER=4

## fine-tuning
export MODEL_PATH=DNABERT/models/${KMER}-new-12w-0
export DATA_PATH=../002_Data_splitting/output_examples/ft_known_${KMER}_${tokens}_${fold}
export OUTPUT_PATH=output_examples/ft_known_output_${KMER}_${tokens}_${fold}_${LR}
export CACHED_PATH=output_examples/cached_dev_${KMER}_${tokens}_${fold}_${LR}


python3 DNABERT/script/run_finetune_20230321.py \
	--model_type dna \
	--tokenizer_name $MODEL_PATH  \
	--model_name_or_path $MODEL_PATH \
	--task_name dnaprom \
	--do_train \
	--do_eval \
	--data_dir $DATA_PATH \
	--cache_dir $CACHED_PATH \
	--max_seq_length ${tokens} \
	--per_gpu_eval_batch_size=64   \
	--per_gpu_train_batch_size=32   \
	--gradient_accumulation_steps 4 \
	--learning_rate ${LR} \
	--num_train_epochs 40 \
	--output_dir $OUTPUT_PATH \
	--evaluate_during_training \
	--logging_steps 1000 \
	--save_steps 10000 \
	--max_steps 10000 \
	--warmup_percent 0.1 \
	--hidden_dropout_prob 0.1 \
	--overwrite_output \
	--weight_decay 0.01 \
	--early_stop 6 \
	--n_process 8 


## predict
export MODEL_PATH=output_examples/ft_known_output_${KMER}_${tokens}_${fold}_${LR}
export DATA_PATH=../002_Data_splitting/output_examples/pred_known_${KMER}_${tokens}_${fold}
export PREDICTION_PATH=output_examples/pred_known_output_${KMER}_${tokens}_${fold}_${LR}
export CACHED_PATH=output_examples/cached_pred_known_${KMER}_${tokens}_${fold}_${LR}

python3 DNABERT/script/run_finetune_20230321.py \
    --model_type dna \
    --tokenizer_name $MODEL_PATH  \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_predict \
    --data_dir $DATA_PATH  \
	--cache_dir $CACHED_PATH \
    --max_seq_length 500 \
    --per_gpu_pred_batch_size=32  \
    --gradient_accumulation_steps 4 \
    --output_dir $MODEL_PATH \
    --predict_dir $PREDICTION_PATH \
    --overwrite_output \
    --n_process 8