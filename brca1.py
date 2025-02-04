import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt
import seaborn as sns

# Загрузка таблиц
df_control = pd.read_csv("MDA_MB_231.txt", sep="\t")  # Контрольные клетки
df_resistant = pd.read_csv("RT_R_MDA_MB_231.txt", sep="\t")  # Радиорезистентные клетки

# Переименование колонок для ясности
df_control.columns = ["Gene", "Control"]
df_resistant.columns = ["Gene", "Resistant"]

# Объединение таблиц по генам
df = pd.merge(df_control, df_resistant, on="Gene")

# Удаление дубликатов и фильтрация генов с нулевой экспрессией в обеих группах
df = df.drop_duplicates(subset=["Gene"])
df = df[(df["Control"] > 0) | (df["Resistant"] > 0)]
print("Очищенные данные:\n", df.head())

# Сохранение очищенной таблицы
df.to_csv("cleaned_gene_expression_data.csv", index=False)
print("Файл cleaned_gene_expression_data.csv сохранен.")

# Расчет Log2FoldChange и p-value
df["Log2FoldChange"] = np.log2((df["Resistant"] + 1) / (df["Control"] + 1))
df["p-value"] = ttest_rel(df["Control"], df["Resistant"])[1]

# Фильтрация значимых генов (p-value < 0.05)
significant_genes = df[df["p-value"] < 0.05]
print("Значимые гены:\n", significant_genes.head())

# Сохранение значимых генов
significant_genes.to_csv("significant_genes.csv", index=False)
print("Файл significant_genes.csv сохранен.")

# Визуализация Volcano Plot
plt.figure(figsize=(8, 6))
plt.scatter(df["Log2FoldChange"], -np.log10(df["p-value"]), c='grey')
plt.scatter(significant_genes["Log2FoldChange"], -np.log10(significant_genes["p-value"]), c='red')
plt.axhline(-np.log10(0.05), color='blue', linestyle='--')
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 p-value")
plt.title("Volcano Plot")
plt.show()

# Heatmap для топ-20 генов
top_genes = significant_genes.nlargest(20, "Log2FoldChange")
plt.figure(figsize=(10, 6))
sns.heatmap(top_genes[["Control", "Resistant"]].set_index(top_genes["Gene"]), cmap="viridis", annot=True)
plt.title("Top 20 Differentially Expressed Genes")
plt.show()
