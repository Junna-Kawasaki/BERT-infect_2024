# ⚠️ Dataset Update

We identified that a subset of the publicly released dataset contained instances that qualify as data leakage.

## Summary of the Issue

Among the 26 virus families evaluated, two families were affected in the time-split datasets.

The leakage involved a limited number of samples in the training data:

- Orthomyxoviridae: 26 out of 8,674 samples (0.3%)

- Rhabdoviridae: 12 out of 127 samples (9.45%)

## Impact on Model Comparison

All models evaluated in the original study were trained and tested using the same dataset configuration.
Therefore, the main conclusion remains unchanged: under identical conditions, BERT-infect demonstrated superior performance compared with existing models.

## Updated Files
- dataset_ver2.xlsx: Corrected dataset with leakage samples removed.

- removed_acc.tsv: Accession numbers of sequences identified as data leakage.

- dataset_ver1.xlsx: Original dataset as used in the published manuscript (archived for reference).
