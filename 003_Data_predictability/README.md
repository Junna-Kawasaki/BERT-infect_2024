# Data predictability 

## Descriptions 
### potential dataset predictability using BLASTn
We used BLASTn to evaluate potential dataset predictability based on the simple hypothesis that a virus has the same infectivity as that with the highest sequence similarity. 


## Scripts  
- **Dataset_predictability_blast.ipynb**: provides an example script for evaluating dataset predictability using BLASTn.  


## Inputs  
- **known_virus_final_dataset_Coronaviridae.txt**  
Here, we provide examples of the family Coronaviridae. The csv and fasta files can be obtained from the Zenodo repository (doi: )


## Outputs
- **blast_prediction_results.csv** 
The column of "subject label" means the predicted infectivity (i.e., animal-infecting (0) or human-infecting (1)).

- **blast_metrics.csv**  
AUROC and PR-AUC scores of BLASTn prediction were listed.
