# Evaluation of model performance

## Descriptions 
### Evaluation of viral infectivity prediction for the past viral genome  
The predictive performances were compared using two metrics, AUROC and PR-AUC. The outputs are associated with Fig 2c in our paper. 

### Evaluation of model detection ability for human-infecting virus using high-throughput sequencing data
We evaluated model predictive performance when using high-throughput sequencing data based on different inputs: (i) 250 bp single-ended high-throughput sequence reads and (ii) variable-length contig sequences (i.e., 500 bp, 1,000 bp, 3,000bp, and 5,000 bp). The predictive performances were compared using two metrics, AUROC and PR-AUC. The output figure is associated with Figure 2d in our paper.

### Evaluation of viral infectivity prediction for the future viral genome
Model predictive performance for the future viral dataset was evaluated under the two scenarios. First, the infectivity of novel viruses was determined based on the threshold representing the highest F1 score in the past viral datasets, and F1, recall and precision metrics were calculated. Second, we investigated the enrichment of human-infecting viruses within the top 20% of predicted probabilities. The output figures are associated with Figures 3b-c, 4b, and 5b-c in our paper.  

## Scripts  
- **Metrics_past_genome_visualization.ipynb**: provides an example script for evaluating model performance for past viral datasets at the genomic sequence level.  

- **Metrics_past_subsequence_visualization.ipynb**: provides an example script for evaluating model performance for short subsequences.  

- **Metrics_future_genome_visualization.ipynb**: provides an example script for evaluating model performance for  viral datasets at the genomic sequence level.  


## Inputs  
- **ft_knwon_4_250_fold_0**: training and evaluating datasets  
- **pred_known_4_250_fold_0**: prediction dataset (past viruses) 
- **pred_known_output_4_250_fold_0**: prediction results for a past dataset
- **pred_unknown_4_250**: prediction dataset (future viruses)
- **pred_unknown_output_4_250_fold_0**: prediction results for a future dataset

## Outputs
The metrics are calculated/visualized when performing each jupyter notebook.