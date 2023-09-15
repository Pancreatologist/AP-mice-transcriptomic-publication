It is the code for the publication titled with 'Circulatory transcriptomic profiling identifies innate-adaptive immunity imbalance in 3 different mouse models of acute pancreatitis' written by Di Wu from in the University of Liverpool. 
The rawdata of microarray is reposed in GEO database, with GEO number depending.
Step:
1.	Use the ‘Library package.R’ to load all the package which we need.
2.	'designwithbatch.txt' is the target file for the readTargets() function in the limma package.
3.	Use the ‘Pre-Process for microarray.R’ to get the gene expression matrix.
4.	‘Preliminary analysis via PCA.R’ is used to make the Figure 3B, Figure 4C and Figure S4 (PCA results).
5.	‘Differentially expressed genes analysis.R’ is used to do the DE analyse to make the Figure 3A, Figure Table S1.
6.	‘Enrichment analyze and visualization (ORA for special genes).R’ is used to do the ORA analyse and make the Figure 3E.
7.	‘GSEA and visualization.R’ is used to make the Figure 3C, Figure 5A and Figure 5B.
8.	'Gene expression heatmap plot' is used to make the heatmap plot for gene expression (Figure 3D, 4E, 5C).
9.	‘Immune inflation.R’ is used to make the gene expression matrix without log transformed and to make the Figure 5D.
10.	'WCGNA.R' is used to do the WCGNA analyse to produce the Figure 4A, 4B and S5.
If you have any comments, please contact me via wudi@liverpool.ac.uk

<!---
Pancreatologist/Pancreatologist is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
