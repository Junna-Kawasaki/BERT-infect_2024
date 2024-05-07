# Dataset curation  
Briefly, our datasets were constructed by (i) collecting virus nucleotide sequences and their metadata from the NCBI Virus database (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/), (ii) linking genomic sequences for each virus strain based on the combination of their metadata, and (iii) excluding redundant sequences with 100% sequence similarity and the identical infectivity labels.

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

## Descriptions  
For segmented RNA viruses:
- Sequences were grouped into viral isolates based on the combination of their metadata.  
- If a single viral isolate contained more sequences than the specified number of viral segments (i.e., Orthomyxoviridae: 8, Rotaviridae: 12, and other segmented RNA viruses: 3), redundancy was eliminated by randomly sampling two sequences for each segment.  

For non-segmented viruses:
- Representative sequences were selected from identical viral sequences with the same host label.  

## Scripts  
- **Dataset_curation_non-segmented_template.ipynb**  
provides an example script for non-segmented viral data curation.  
- **Dataset_curation_segmented_template.ipynb**   
provides an example script for segmented viral data curation.  


## Inputs  
- **${virus_family}.csv**  
Viral metadata obtained from the NCBI Virus Database. Here, we provide examples of the families Coronaviridae and Hantaviridae.


## Outputs
- **${virus_family}.curated.csv**  
The complete files for 26 viral families are published as the supplementary table 2.