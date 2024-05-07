# fine-tuning BERT models  

## Descriptions 
We constructed new predictive models for viral infectivity, namely BERT-infect, based on the LLMs. We used the (i) DNABERT model pre-trained on human whole genome and (ii) ViBE model pre-trained on the viral genome sequences in the NCBI RefSeq database. 


## Scripts  
- **DNABERT_fintuning.sh**: provides an example script for fine-tuning the DNABERT model.  
- **ViBE_fintuning.sh**: provides an example script for fine-tuning the ViBE model.  

- **run_finetune_20230321.py** is a modified version of published script in the following respects:
    - To address class imbalance, we modified the loss function to a weighted binary cross-entropy. Please replace “DNABERT/examples/run_finetune.py” with this script.

## Environment setup
Please check the DNABERT git-hub repository (https://github.com/jerryji1993/DNABERT) and ViBE git-hub repository (https://github.com/DMnBI/ViBE).

## Inputs
The example of inputs are placed in *"002_Data_splitting/output_examples"*.

## Outputs
The BERT-infect models fine-tuned for each virus family are published in the Zendo repository:
- BERT-infect models for non-segmented DNA viruses (doi: 10.5281/zenodo.11103056)
- BERT-infect models for non-segmented RNA viruses (doi: 10.5281/zenodo.11103079)
- BERT-infect models for segmented RNA viruses (doi: 10.5281/zenodo.11103091)