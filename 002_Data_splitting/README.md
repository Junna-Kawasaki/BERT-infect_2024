# Data splitting 
To investigate model generalization capability, we divide viral data into two datasets: past and future viral datasets.

The dataset information is published as the supplementary table 2.  

## Descriptions 
### Past and future virus datasets 
The collected data were divided into two datasets: (i) a past virus dataset for model training and (ii) a future virus dataset for evaluating model generalization capability. Data with a sequence collection date before December 31, 2017, were classified into the past virus dataset, while the subsequent data were classified into the future virus dataset.   

### five-fold stratified cross-validation
The past virus datasets were divided for five-fold stratified cross-validation. The fold datasets were prepared to maintain the same proportions of infectivity labels and viral genera as those in the overall data.

## Prepare sequence data
The curated datasets (i.e., excel file for metadata and fasta/genbank file for nucleotide sequences) are available in the Zenodo repository (doi: 10.5281/zenodo.11102793). Please download the fasta and metadata files from Zenodo and prepare as inputs.  

```
# download fasta file
wget \
    https://zenodo.org/records/11102793/files/Coronaviridae.tar.xz \
    tar Jxfv Coronaviridae.tar.xz 

cat Coronaviridae/future.fasta \
    Coronaviridae/past.fasta \
    > input_examples/Coronaviridae.curated.fasta
```

## Scripts  
- **Data_splitting_5CV_template.ipynb**: provides an example script for data splitting for five-fold stratified cross-validation.  


## Inputs  
- **${virus_family}.curated.csv**  
- **${virus_family}.curated.fasta**  

## Outputs
- **/ft_known_4_250_${fold}/train.tsv**
- **/ft_known_4_250_${fold}/dev.tsv**
- **/pred_known_4_250_${fold}/dev.tsv**

For each fold dataset, three files for training, evaluating, and testing datasets were created.
