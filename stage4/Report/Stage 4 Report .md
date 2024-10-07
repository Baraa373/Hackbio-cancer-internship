<!--StartFragment-->

**Exploring Gene Expression Patterns in Gliomas: Clustering Approaches for IDH Status Determination**

**Authors:**  Braa Elwaleed Mohamed (@Baraa),

Nooran\_ Alharty (@ Nooran\_ Alharty), Henry Momanyi (@Henry), Taha (@Taha),Olabode Kaosara Omowunmi (@TheResearchNerd), Olajumoke Oladosu (@ jumblosu).

**Introduction**

Glioma is a kind of brain tumor that arises from glial cells. They can be categorized based on their mutations into different subtypes, ranging from low-grade gliomas (LGG) to high-grade gliomas such as glioblastomas. IDH1 and IDH2 are the most significantly linked to a better prognosis than IDH-wildtype gliomas. Methylation profiling has been explained to differentiate between IDH-mutant and IDH-wildtype gliomas.

**Description** 

The TCGA LGG dataset of 516 samples focused on IDH mutation status was processed using TCGAbiolinks.After removing samples with missing IDH data, gene expression was log2-transformed and filtered for clustering and biomarker discovery.

**Methodology**

The data collection and analysis process is Represented in the workflow diagram below. 

<!--StartFragment-->

**Figure 1: Workflow Overview**

**Results**

The KNN classifier, trained on selected features with cross-validation, achieved 79% accuracy in testing. The model successfully classified most samples as either IDH mutant or wild type.

<!--StartFragment-->

**figure 2: KNN Confusion Matrix**

The heatmap displays cluster means of gene expression data from TCGA samples. It groups genes and samples based on expression levels, using color intensity to indicate higher or lower expression. This clustering helps visualize relationships and patterns in gene expression across different samples.

<!--StartFragment-->

Figure 3: Heatmap of cluster means

<!--EndFragment-->

\


<!--EndFragment-->

<!--EndFragment-->

<!--EndFragment-->

<!--StartFragment-->

Figure 5: upregulation gene

We have upregulated gene pathways. The TGF-beta signaling pathway has an FDR of 3.9E-03 and a fold enrichment of 5.3 (93 genes). The Glutamatergic synapse has an FDR of 1.3E-02 and a fold enrichment of 4.3 (114 genes). The Neuroactive ligand-receptor interaction shows an FDR of 3.9E-05 and a fold enrichment of 3.6 (350 genes). The CAMP signaling pathway has an FDR of 3.6E-03 and a fold enrichment of 3.5 (219 genes).

<!--StartFragment-->

Figure 4: Downregulation gene

Downregulation gene pathway enrichment results using ShinGo. The ECM-Receptor Interaction pathway has an FDR of 5.9E-13 and a fold enrichment 6.3. The Primary Immunodeficiency pathway shows an FDR of 2.6E-06 with a fold enrichment of 4.3. The Basal Cell Carcinoma pathway has an FDR of 1.5E-05 and a fold enrichment of 3.4.

**Discussion****** 

**DEA** analysis with a Prediction Forest model showed 79% accuracy in IDH detection, indicating higher expression in mutant gliomas (frequency 120) compared to wildtype (30%). 

After filtering to 100 genes, the heat map revealed \*\*two significant clusters that effectively separate the samples. While more clusters could refine the analysis, the current two clusters demonstrate clear segregation.

**Conclusion** 

More data collection is needed to improve IDH classification and treatment strategies in glioma research.

**References** 

Kellerw, S. A., Gupta, A., & Hsu, F. H. (2015). Mutant IDH1 and IDH2 in glioma: A review of their roles in the development and progression of brain tumors. \*Cell\*, \*162\*(4), 901-914. https\://doi.org/10.1016/j.cell.2015.08.023

\


<!--EndFragment-->

<!--EndFragment-->
