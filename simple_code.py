#!/usr/bin/env python
# coding: utf-8

# You are going to complete the following tasks to reinforce your python skills in bioinformatics.
# 
# - Write a python function for translating DNA to protein

# In[1]:


from Bio.Seq import Seq

def dna_to_protein(DNAseq):
    #Clean the input
    DNAseq_clean = DNAseq.replace('\n','').replace(' ','').upper().strip()
    
    #Validate the sequence
    valid_bases = {'A', 'T', 'C', 'G'}
    if not DNAseq_clean:
        print("Error: DNA sequence is empty.")
        return None
    elif any(base not in valid_bases for base in DNAseq_clean):
        print("Error: DNA sequence contains invalid characters. Only A, T, C, G are allowed.")
        return None
    else:
        # Step 3: Create a Seq object and translate
        dna_seq = Seq(DNAseq_clean)
        protein_seq = dna_seq.translate(to_stop=True)
        return str(protein_seq)


# - Write a python function for calculating the hamming distance between your slack username (e.g josoga) and twitter/X (joseph) handle (synthesize one if you don’t have one). Feel free to pad it with extra words if they are not of the same length.

# In[2]:


def hamming_distance(str1, str2):
    #Pad the shorter string with spaces
    len1 = len(str1)
    len2 = len(str2)
    if len1 < len2:
        str1 = str1 + " " * (len2 - len1)
    elif len2 < len1:
        str2 = str2 + " " * (len1 - len2)
    
    #Count the differences
    distance = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            distance += 1
    print(str1)
    print(str2)
    return distance

slack_username = "debangana"
twitter_handle = "debanganabala"

distance = hamming_distance(slack_username, twitter_handle)
print("Hamming distance:", distance)


# # PART A GENE EXPRESSION ANALYSIS ##############

# In[3]:


#Required libraries 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np


# 1. a. Heatmap
# - Use the normalized gene expression dataset to plot a clustered heatmap of the top differentially expressed genes between HBR and UHR samples.
# 
# - Label both genes and samples.
# 
# - Use a color gradient (e.g., Blues) to indicate expression levels.

# In[4]:


# Load normalized gene expression data (genes as rows, samples as columns)
df = pd.read_csv("hbr_uhr_top_deg_normalized_counts.csv", index_col=0)

# Create clustered heatmap with normalization and blue color gradient
g = sns.clustermap(
    df, cmap="Blues", standard_scale=1, figsize=(10, 8),
)

# Add centered title and display
g.fig.suptitle("Clustered Heatmap of Top Differentially Expressed Genes between HBR and UHR samples", y=1.02)
plt.show()


# 2. b. Volcano Plot
# 
# - Plot log2FoldChange vs log10(Padj) from the DEG results.
# - Color points by significance:
# - Upregulated: green
# - Downregulated: orange
# - Not significant: grey
# - Add dashed vertical lines at log2FoldChange = ±1

# In[5]:


# Load DEG results
hbr_uhr_deg_chr22 = pd.read_csv("hbr_uhr_deg_chr22_with_significance.csv")

#creating 
plot1 = sns.scatterplot(
    data=hbr_uhr_deg_chr22,
    x="log2FoldChange",
    y="-log10PAdj",
    hue="significance"
)


# Add vertical dashed lines at ±1
plt.axvline(x=1, color="gray", linestyle="--", linewidth=1)
plt.axvline(x=-1, color="gray", linestyle="--", linewidth=1)

# Add horizontal dashed line at -log10(0.05) ≈ 1.3
plt.axhline(y=1.3, color="gray", linestyle="--", linewidth=1)

# Adding title to the Volcano Plot
plt.title("Volcano Plot (chr22)")
sns.despine()


# # Part B – Breast Cancer Data Exploration

# 1.c. Scatter Plot (radius vs texture)
# - Plot texture_mean vs radius_mean and color points by diagnosis (M = malignant, B = benign).

# In[6]:


data = pd.read_csv("data-3.csv")

plt.figure()
sns.scatterplot(data=data, x='radius_mean', y='texture_mean', hue='diagnosis', palette={'M':'red', 'B':'blue'})
plt.xlabel("radius_mean")
plt.ylabel("texture_mean")
plt.show()


# 2.d. Correlation Heatmap
# - Compute the correlation matrix of six key features:
# 
# radius_mean, texture_mean, perimeter_mean, area_mean, smoothness_mean, compactness_mean.
# 
# - Plot as a heatmap with correlation values annotated.

# In[7]:


features = ['radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean']
corr = data[features].corr()
plt.figure()
sns.heatmap(corr, annot=True, cmap='coolwarm')
plt.show()


# 3. e. Scatter Plot (smoothness vs compactness)
# - Plot compactness_mean vs smoothness_mean colored by diagnosis.
# - Include gridlines and clear axis labels.

# In[8]:


plt.figure()
sns.scatterplot(data=data, x='smoothness_mean', y='compactness_mean', hue='diagnosis', palette={'M':'red', 'B':'blue'})
plt.xlabel("smoothness_mean")
plt.ylabel("compactness_mean")
plt.grid(True)
plt.show()


# 4.f. Density Plot (area distribution)
# - Plot kernel density estimates (KDE) of area_mean for both M and B diagnoses on the same axis.
# - Add legend and labeled axes.

# In[9]:


m = data[data['diagnosis']=='M']['area_mean']
b = data[data['diagnosis']=='B']['area_mean']
x = np.linspace(min(m.min(), b.min()), max(m.max(), b.max()), 500)
plt.plot(x, gaussian_kde(m)(x), label='M')
plt.plot(x, gaussian_kde(b)(x), label='B')
plt.xlabel("area_mean")
plt.ylabel("Density")
plt.legend()
plt.show()

