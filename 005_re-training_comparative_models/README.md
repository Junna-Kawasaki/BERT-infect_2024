# Re-training comparative models

## Descriptions 
Four existing models (humanVirusFinder, DeePac_vir, VIDHOP, and zoonotic_rank) were re-trained with our newly constructed datasets. The hyperparameters for the DeePac_vir and VIDHOP model re-training are listed in Supplementary Table 6. The humanVirusFinder and zoonotic_rank models were re-trained with the same parameters as in the previous studies.  

## Example scripts  
- **DeePac_vir_train_pred.sh**
- **VIDHOP_train_pred.sh**
- **humanVirusFinder_train_pred.sh**
- **zoonotic_rank_train_pred.sh**

## Environment setup  
Before running example scripts, please setup your environment following the git-hub repositories.

| model | git-hub |  
| --- | --- |
| DeePac_vir| https://gitlab.com/dacs-hpi/deepac-vir |  
| humanVirusFinder | https://github.com/Fzhang1992/human-infecting_virus_finder |  
| VIDHOP | https://github.com/flomock/vidhop |
| zoonotic_rank | https://github.com/Nardus/zoonotic_rank |
