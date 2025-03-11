# RNA-seq Differential Expression Analysis

## Overview
This project performs differential expression analysis on RNA-seq data, comparing gene expression levels between control and radioresistant conditions. It includes:
- Data loading and preprocessing
- Calculation of Log2FoldChange and p-values
- Multiple testing correction using the Benjamini-Hochberg method
- DESeq2 integration for differential expression analysis
- Visualization via Volcano Plot and Heatmap

## Requirements
### Python Packages
Make sure you have the following Python packages installed:
```sh
pip install pandas numpy scipy matplotlib seaborn statsmodels tkinter rpy2
```

### R Dependencies
This project integrates with DESeq2 in R. Ensure you have R installed along with the required packages:
```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Usage
1. Run the script:
```sh
python script.py
```
2. Select the control and radioresistant dataset files when prompted.
3. The script will process the data, compute differentially expressed genes, and generate visualizations.
4. Results will be saved as:
   - `cleaned_gene_expression_data.csv` – Processed gene expression data
   - `significant_genes.csv` – Filtered significant genes
   - Volcano and heatmap plots displayed

## File Format
The input files should be tab-separated text files with two columns:
```
Gene    Expression
GENE1   50
GENE2   75
...
```

## Output Files
- `cleaned_gene_expression_data.csv`: Contains gene expression values, Log2FoldChange, and adjusted p-values.
- `significant_genes.csv`: Filtered genes with an adjusted p-value < 0.05.

![alt text](https://github.com/AkimzhanRakhimov/GeneExpressionAnalysis/blob/main/2025-03-11_14-29-27.png)

## Visualization
- **Volcano Plot**: Highlights significantly differentially expressed genes.
- **Heatmap**: Shows expression levels of the top 20 differentially expressed genes.

## License
This project is open-source and available under the MIT License.

