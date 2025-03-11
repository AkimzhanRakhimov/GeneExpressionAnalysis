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
# Функции для анализа данных RNA-seq
# ===============================

def load_and_prepare_data():
    # Открытие окна выбора файлов
    Tk().withdraw()  # Скрытие главного окна Tk
    print("Выберите файл с контрольными данными")
    control_file = filedialog.askopenfilename(title="Выберите файл с контрольными данными")
    print("Выберите файл с радиорезистентными данными")
    resistant_file = filedialog.askopenfilename(title="Выберите файл с радиорезистентными данными")

    if not control_file or not resistant_file:
        print("Ошибка: не выбраны файлы для анализа.")
        return None

    # Загрузка данных
    df_control = pd.read_csv(control_file, sep="\t")
    df_resistant = pd.read_csv(resistant_file, sep="\t")

    # Переименование колонок
    df_control.columns = ["Gene", "Control"]
    df_resistant.columns = ["Gene", "Resistant"]

    # Объединение таблиц и фильтрация данных
    df = pd.merge(df_control, df_resistant, on="Gene").drop_duplicates(subset=["Gene"])
    df = df[(df["Control"] > 0) | (df["Resistant"] > 0)]

    print("Данные успешно загружены и подготовлены.")
    return df

# Функция расчета Log2FoldChange и p-value с коррекцией методом Бенджамини-Хохберга
def compute_differential_expression(df):
    df["Log2FoldChange"] = np.log2((df["Resistant"] + 1) / (df["Control"] + 1))
    df["p-value"] = ttest_rel(df["Control"], df["Resistant"])[1]

    # Коррекция p-value методом Бенджамини-Хохберга
    df["adjusted p-value"] = multipletests(df["p-value"], method="fdr_bh")[1]
    return df

# Функция фильтрации значимых генов
def filter_significant_genes(df, p_threshold=0.05):
    significant_genes = df[df["adjusted p-value"] < p_threshold]
    print(f"Найдено {len(significant_genes)} значимых генов (adjusted p-value < {p_threshold}).")
    return significant_genes

# Использование DESeq2 для анализа дифференциальной экспрессии через R
def compute_deseq2_differential_expression(df):
    try:
        deseq2 = importr("DESeq2")
        base = importr("base")

        # Создание объекта DataFrame для R
        ro.globalenv["count_data"] = pandas2ri.py2rpy(df.set_index("Gene"))
        design_matrix = ro.r("data.frame(condition=factor(c(rep('Control', nrow(count_data)/2), rep('Resistant', nrow(count_data)/2))))")

        # Анализ через DESeq2
        dds = deseq2.DESeqDataSetFromMatrix(countData=ro.globalenv["count_data"], colData=design_matrix, design=ro.Formula("~ condition"))
        dds = deseq2.DESeq(dds)
        results = deseq2.results(dds)

        # Преобразование результата обратно в DataFrame
        res_df = pandas2ri.rpy2py(base.as_data_frame(results))
        res_df.columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "p-value", "adjusted p-value"]

        print("DESeq2 анализ завершен.")
        return res_df
    except Exception as e:
        print("Ошибка при использовании DESeq2:", e)
        return df

# Функция сохранения данных в CSV
def save_data(df, filename):
    df.to_csv(filename, index=False)
    print(f"Данные сохранены в файл {filename}.")

# Функция визуализации Volcano Plot
def plot_volcano(df, significant_genes):
    plt.figure(figsize=(8, 6))
    plt.scatter(df["Log2FoldChange"], -np.log10(df["p-value"]), c='grey')
    plt.scatter(significant_genes["Log2FoldChange"], -np.log10(significant_genes["p-value"]), c='red')
    plt.axhline(-np.log10(0.05), color='blue', linestyle='--')
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 p-value")
    plt.title("Volcano Plot")
    plt.show()

# Функция визуализации Heatmap для топ-N генов
def plot_heatmap(significant_genes, top_n=20):
    top_genes = significant_genes.nlargest(top_n, "Log2FoldChange")
    plt.figure(figsize=(10, 6))
    sns.heatmap(top_genes[["Control", "Resistant"]].set_index(top_genes["Gene"]), cmap="viridis", annot=True)
    plt.title(f"Top {top_n} Differentially Expressed Genes")
    plt.show()

# ===============================
# Основной анализ
# ===============================

def main():
    # Загрузка и обработка данных
    df = load_and_prepare_data()
    if df is None:
        return

    # Расчет дифференциальной экспрессии
    df = compute_differential_expression(df)

    # Сохранение очищенных данных
    save_data(df, "cleaned_gene_expression_data.csv")

    # Фильтрация значимых генов
    significant_genes = filter_significant_genes(df)

    # Сохранение значимых генов
    save_data(significant_genes, "significant_genes.csv")

    # Визуализации
    plot_volcano(df, significant_genes)
    plot_heatmap(significant_genes)

if __name__ == "__main__":
    main()