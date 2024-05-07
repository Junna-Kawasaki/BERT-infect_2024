# Phylogenetic analysis 

## Descriptions 
Multiple sequence alignments for each viral family were constructed by adding sequences to the reference alignment provided by the International Committee on Taxonomy of Viruses resources.  
Phylogenetic trees were constructed by the maximum likelihood method using IQ-TREE (version 2.1.4 beta). Substitution models were selected based on the Bayesian information criterion provided by ModelFinder.  

```
wget https://zenodo.org/records/11102793/files/Flaviviridae.tar.xz
tar Jxfv Flaviviridae.tar.xz 

cat Flaviviridae/future.fasta \
    Flaviviridae/past.fasta \
    > input_examples/Flaviviridae.curated.fasta
```


## Scripts  
- **Phylogenetic_tree_aa_corona.ipynb**: provides an example script for making tree of Coronaviridae.  

- **Phylogenetic_tree_visualization_corona.ipynb**: provides an example script for visualizing tree of Coronaviridae. The output figure corresponds to Fig 4a.

- **Phylogenetic_tree_nucl_flavi.ipynb**: provides an example script for making tree of Flaviviridae.  

- **Phylogenetic_tree_visualization_flavi.ipynb**: provides an example script for visualizing tree of Flaviviridae. The output figure corresponds to Fig 5a.

## Inputs  
 - *Coronaviridae* and *Flaviviridae* folders contains tree files. The original data for phylogenetic analysis were uploaded in the Zenodo repository (doi: )

## Outputs
The visualization is conducted by each jupyter notebook.