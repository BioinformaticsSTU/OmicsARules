# OmicsARules
a R package for integration of multi-omics datasets via association rules mining


## Introduction
OmicsARules is a tool for the analysis of multi-omics high-throughput data based on the use of association rules. OmicsARules supports to identify recurrent and associated patterns, and provides a new dimension for exploring single or multiple omics data across sequencing platforms or across samples. Besides, a new rule-interestingness measure Lamda3 was Embeded, it can be used to evaluate the association rules and identify biologically significant patterns. Association rule mining and visualizing were implemented in R environment using package arules and ggplot2.

## Installation

Dependencies should be installed firstly.
```
install.packages("arules");
install.packages("tcltk");
install.packages("ggplot2");
install.packages("igraph");
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager");
BiocManager::install("GOSim");
```

After requirement is satisfied, you can install omicsArules by running command below.
```
install.packages("omicsArules_0.0.1.tar.gz");
```

## What is lamda3?

Lamda3 is a new measure to indicates the significance and interest of rules. This measure not only use binary data, but also continuous variables to overcome the shortcoming of loss of information that other association rules mining algorithms(ARM) used to have.

<img src="img/LAMDA.png"  width="800" height="900" />

## Data discretization

The uploading data is preferable to be a binary dataset containing interested genes. Alternatively, if the input dataset contains continuous variable, users should firstly identify the interesting genes according to their own measures. For example, for mRNA profiling data, the genes can be selected and sorted by P values from differential expression analysis. OmicsARules provides four cutoff values to discretize the continuous values into binary matrix, namely mean, median, P25 (the upper quartile), P75 (the lower quartile). OmicsARules calculates one of these cutoff values (according to the user's choice) in each column, and if the values in each gene of a particular sample larger than the cutoff value, this value would be transformed into 1, otherwise, 0 is used. 

## Result example

In this example, association rules analysis was performed on ESCA mRNA expression (support=0.4 and confidence=0.8). The dataset contains 184 patients and the top-2000 DEGs. Interestingly, some well-known relationships between genes were observed. For instance, a particular rule {TAP1}==>{PSMB9} was identified with support 0.413, confidence 0.883, lift 1.711 and lamda3 21.445. From biological viewpoint, this rule means co-occurrence of TAP1 and PSMB9 dysregulation in mRNA expression happened on more than 40% ESCA patients [the actual frequency: supp(TAP1 ∪ PSMB9 )=41.3%]. In view of confidence, when TAP1 was dysregulated, the possibility of simultaneously altered PSMB9 expression was 88.3% . As for another measure lift, compared to possibility of random events, that is, dysregulation of TAP1 was independent of CCNB1 , their co-dysregulation was 1.711 times more frequent.

```
 lhs rhs support confidence lift lamda3 
1 {ZNF329} {ZNF793} 0.347 0.82 1.65 11.49  
2 {TAP1} {PSMB9} 0.413 0.883 1.711 21.445  
3 {PSMB9} {TAP1} 0.413 0.8 1.711 14.723  
4 {PSME1} {PSMB9} 0.402 0.831 1.61 9.263  
5 {ZNF345} {ZNF790} 0.385 0.855 1.873 41.334  
```

## Mining association rule conserved across multiple omics datasets

Regarding multi-omics datasets, for example, both mRNA expression and DNA methylation datasets retrieved from the same cohort, each dataset was separately subjected to the preprocess step such as differential expression analysis and discretization. Then these two binary datasets were combined according the sample IDs, and then subjected to OmicsARules association mining. In order to discriminate sources of genes, suffix '.1' or '.2' was added behind the gene symbol, while the former indicated genes were present in the mRNA dataset; and the latter indicated genes were from the DNA methylation. Finally, rule-interestingness measures, namely Lift and Lamda3, as well as Wcs (weighted condensed support), were calculated to rank the rules.

```
lhs	rhs	support	confidence	lift	lamda3
1	{CDKN3.1,A1BG.2}	{CDC20.1}	0.396	0.811	1.448	0
2	{CDKN3.1,A1CF.2}	{CDC20.1}	0.396	0.811	1.448	0
3	{NEK2.1,AADACL3.2}	{UBE2T.1}	0.423	0.804	1.541	0
4	{UBE2T.1,AADACL2.2}	{NEK2.1}	0.423	0.812	1.541	0
5	{KIF2C.1,AADACL3.2}	{CDC20.1}	0.456	0.848	1.515	0
6	{CDC20.1,AADACL3.2}	{KIF2C.1}	0.456	0.815	1.515	0
7	{CDKN3.1,AADAT.2}	{CDC20.1}	0.396	0.811	1.448	0

```
An example of the common rules.

## Group based demonstration 

To visualize the grouped matrix, a balloon plot was created with antecedent groups as columns (LHS) and consequents as rows (RHS). The idea is that genes on the left side of several rules, which are statistically dependent on the same gene on the right side, are supposed to be similar and thus can be grouped together. We start with the set of association rules: R = { a1 , c1 , m1 , . . . ai , ci , mi , . . . an , cn , mn }, where ai is the gene or gene set on the LHS, ci is the gene on the RHS and mi is the selected interest measure (default: lift) for the i-th rule for i = 1, . . . , n. Consequently, a L × K matrix M in R with one column for each unique antecedent and one row for each unique consequent was identified and created. The color of balloons represent the lift of rules and the size of balloons represent the support of rules. For instance, one of column names, '3 rules: {HS3ST4, KCNMB1}', indicated that this particular rule with either one or two genes on the LHS was appeared in 3 rules, with three different genes on the RHS (namely SYNPO2, LMOD1 and MYH11). 

<img src="img/group.png"  width="700" height="450" />
Group plot of association rules

## Graph based demonstration  
This method visualizes association rules using vertices typically represent genes or gene sets and edges indicate their relationship in rules.

<img src="img/graph.png"  width="700" height="450" />
Graph plot of association rules



----------------------------
# License
----------------------------

OmicsARules is released under the MIT license.

