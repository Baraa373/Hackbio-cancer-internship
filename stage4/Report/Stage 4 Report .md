<!--StartFragment-->

**Exploring Gene Expression Patterns in Gliomas: Clustering Approaches for IDH Status Determination**

**Authors:**  Braa Elwaleed Mohamed (@Baraa), Nooran\_ Alharty (@ Nooran\_ Alharty), Henry Momanyi (@Henry), Taha (@Taha), Olabode Kaosara Omowunmi (@TheResearchNerd), Olajumoke Oladosu (@ jumblosu).

**Introduction**

Glioma is a kind of brain tumor that arises from glial cells. They can be categorized based on their mutations into different subtypes, ranging from low-grade gliomas (LGG) to high-grade gliomas such as glioblastomas. IDH1 and IDH2 are the most significantly linked to a better prognosis than IDH-wildtype gliomas. Methylation profiling has been explained to differentiate between IDH-mutant and IDH-wildtype gliomas.

**Description** 

The TCGA LGG dataset of 516 samples focused on IDH mutation status was processed using TCGAbiolinks.After removing samples with missing IDH data, gene expression was log2-transformed and filtered for clustering and biomarker discovery.

**Methodology**

The data collection and analysis process is Represented in the workflow diagram below. 

![flowchart](https://github.com/user-attachments/assets/bc95b38a-c30f-4bf2-a501-69f1b562aa5e)


**Figure 1: Workflow Overview**

**Results**

The KNN classifier, trained on selected features with cross-validation, achieved 79% accuracy in testing. The model successfully classified most samples as either IDH mutant or wild type.


![Rplot01](https://github.com/user-attachments/assets/de8894df-cfcf-4fda-b3b8-1f85ba9c7efd)

**Figure 2: KNN Confusion Matrix**

The heatmap displays cluster means of gene expression data from TCGA samples. It groups genes and samples based on expression levels, using color intensity to indicate higher or lower expression. This clustering helps visualize relationships and patterns in gene expression across different samples.


![Heatmap](https://github.com/user-attachments/assets/85c56d19-4710-46dc-b812-1e98ea145fa3)


Clustering revealed gene groups with a similar expression. This suggests they work together and may indicate important markers for IDH status and tumor behavior.

Figure 3: Heatmap of cluster means


![upregulated](https://github.com/user-attachments/assets/3d15e973-424f-4438-8823-b0925d7452da)

Figure 4: upregulation gene

We have upregulated gene pathways. The TGF-beta signaling pathway has an FDR of 3.9E-03 and a fold enrichment of 5.3 (93 genes). The Glutamatergic synapse has an FDR of 1.3E-02 and a fold enrichment of 4.3 (114 genes). The Neuroactive ligand-receptor interaction shows an FDR of 3.9E-05 and a fold enrichment of 3.6 (350 genes). The CAMP signaling pathway has an FDR of 3.6E-03 and a fold enrichment of 3.5 (219 genes).


![downregulated](https://github.com/user-attachments/assets/47b931d1-f9f6-40c5-af79-ecf03a4fe6a8)

Figure 4: Downregulation gene

Downregulation gene pathway enrichment results using ShinGo. The ECM-Receptor Interaction pathway has an FDR of 5.9E-13 and a fold enrichment 6.3. The Primary Immunodeficiency pathway shows an FDR of 2.6E-06 with a fold enrichment of 4.3. The Basal Cell Carcinoma pathway has an FDR of 1.5E-05 and a fold enrichment of 3.4.

**Discussion**

**DEA** analysis with a Prediction KNN model showed 79% accuracy in IDH detection, indicating higher expression in mutant gliomas (frequency 120) compared to wildtype (30%). 

After filtering to 100 genes, the heat map revealed two significant clusters that effectively separated the samples. While more clusters could refine the analysis, the current two clusters demonstrate clear segregation.

**Conclusion** 

IDH mutations can be identified using methylation profiling and gene expression data. The KNN classifier reached 79% accuracy. More data is needed to improve IDH classification and treatment.

**References** 

Michele Ceccarelli1,2,24 ∙ Floris P. Barthel3,4,24 ∙ Tathiane M. Malta5,6,24∙ … ∙ Houtan Noushmehr5,6,25 houtan@usp.br ∙ Antonio Iavarone10,25 ai2102@columbia.edu ∙ Roel G.W. Verhaak3,25. Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Volume 164, Issue 3p550-563January 28, 2016.https://www.cell.com/cell/fulltext/S0092-8674(15)01692-X


