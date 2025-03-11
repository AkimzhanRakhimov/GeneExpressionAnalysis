import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt
import seaborn as sns
from tkinter import Tk, filedialog
from statsmodels.stats.multitest import multipletests
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

pandas2ri.activate()

# ===============================
# RNA-seq Data Analysis Functions
# ===============================

def load_and_prepare_data():
    # Open file selection dialogs
    Tk().withdraw()  # Hide the main Tk window
    print("Select the control data file")
    control_file = filedialog.askopenfilename(title="Select the control data file")
    print("Select the radioresistant data file")
    resistant_file = filedialog.askopenfilename(title="Select the radioresistant data file")

    if not control_file or not resistant_file:
        print("Error: No files selected for analysis.")
        return None

    # Load data
    df_control = pd.read_csv(control_file, sep="\t")
    df_resistant = pd.read_csv(resistant_file, sep="\t")

    # Rename columns
    df_control.columns = ["Gene", "Control"]
    df_resistant.columns = ["Gene", "Resistant"]

    # Merge tables and filter data
    df = pd.merge(df_control, df_resistant, on="Gene").drop_duplicates(subset=["Gene"])
    df = df[(df["Control"] > 0) | (df["Resistant"] > 0)]

    print("Data successfully loaded and prepared.")
    return df

# Function to compute Log2FoldChange and p-value with Benjamini-Hochberg correction
def compute_differential_expression(df):
    df["Log2FoldChange"] = np.log2((df["Resistant"] + 1) / (df["Control"] + 1))
    df["p-value"] = ttest_rel(df["Control"], df["Resistant"])[1]

    # Adjust p-values using Benjamini-Hochberg method
    df["adjusted p-value"] = multipletests(df["p-value"], method="fdr_bh")[1]
    return df

# Function to filter significant genes
def filter_significant_genes(df, p_threshold=0.05):
    significant_genes = df[df["adjusted p-value"] < p_threshold]
    print(f"Found {len(significant_genes)} significant genes (adjusted p-value < {p_threshold}).")
    return significant_genes

# Using DESeq2 for differential expression analysis via R
def compute_deseq2_differential_expression(df):
    try:
        deseq2 = importr("DESeq2")
        base = importr("base")

        # Create an R DataFrame
        ro.globalenv["count_data"] = pandas2ri.py2rpy(df.set_index("Gene"))
        design_matrix = ro.r("data.frame(condition=factor(c(rep('Control', nrow(count_data)/2), rep('Resistant', nrow(count_data)/2))))")

        # Perform DESeq2 analysis
        dds = deseq2.DESeqDataSetFromMatrix(countData=ro.globalenv["count_data"], colData=design_matrix, design=ro.Formula("~ condition"))
        dds = deseq2.DESeq(dds)
        results = deseq2.results(dds)

        # Convert results back to DataFrame
        res_df = pandas2ri.rpy2py(base.as_data_frame(results))
        res_df.columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "p-value", "adjusted p-value"]

        print("DESeq2 analysis completed.")
        return res_df
    except Exception as e:
        print("Error using DESeq2:", e)
        return df

# Function to save data to CSV
def save_data(df, filename):
    df.to_csv(filename, index=False)
    print(f"Data saved to {filename}.")

# Function to visualize Volcano Plot
def plot_volcano(df, significant_genes):
    plt.figure(figsize=(8, 6))
    plt.scatter(df["Log2FoldChange"], -np.log10(df["p-value"]), c='grey')
    plt.scatter(significant_genes["Log2FoldChange"], -np.log10(significant_genes["p-value"]), c='red')
    plt.axhline(-np.log10(0.05), color='blue', linestyle='--')
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 p-value")
    plt.title("Volcano Plot")
    plt.show()

# Function to visualize Heatmap for top-N genes
def plot_heatmap(significant_genes, top_n=20):
    top_genes = significant_genes.nlargest(top_n, "Log2FoldChange")
    plt.figure(figsize=(10, 6))
    sns.heatmap(top_genes[["Control", "Resistant"]].set_index(top_genes["Gene"]), cmap="viridis", annot=True)
    plt.title(f"Top {top_n} Differentially Expressed Genes")
    plt.show()

# ===============================
# Main Analysis
# ===============================

def main():
    # Load and process data
    df = load_and_prepare_data()
    if df is None:
        return

    # Compute differential expression
    df = compute_differential_expression(df)

    # Save cleaned data
    save_data(df, "cleaned_gene_expression_data.csv")

    # Filter significant genes
    significant_genes = filter_significant_genes(df)

    # Save significant genes
    save_data(significant_genes, "significant_genes.csv")

    # Visualizations
    plot_volcano(df, significant_genes)
    plot_heatmap(significant_genes)

if __name__ == "__main__":
    main()
