[![Creative Commons License](http://i.creativecommons.org/l/by/4.0/80x15.png)](https://github.com/vntasis/SLE_spatial_gene_expression/blob/master/LICENSE)

# Extensive fragmentation and re-distribution of gene co-expression patterns underlie the progression of Systemic Lupus Erythematosus


Repository for the code and the data produced during the analysis described in the article "Extensive fragmentation and re-distribution of gene co-expression patterns underlie the progression of Systemic Lupus Erythematosus"

[Description of the data](#description) | [Description of the code](#code) | [Cite](#cite)

## Description of the data<a name="description"></a>
In this work, we have ananlyzed an RNAseq dataset (Panoussis et al., 2019)
of a large cohort of Systemic Lupus Erythematosus (SLE) patients and
a control group of healthy individuals. We studied gene expression 
from an "architectural" point of view, by measuring gene expression
correlation in a chromosome-wise manner. We statistically evaluated
these correlation values by implementing a robust permutation test 
(x1000). Based on that, we detected 'Domains of Co-ordinated 
expression' (DCEs), as genomic regions with high correlation values
in the interior part and lower correlation with their respective 
neighbourhood.

We performed this analysis for the healthy group and three SLE groups
of increasing disease activity/severity. The results indicate an
extensive reorganization of co-expression patterns in
patients compared to healthy individuals. We further classified these
changes into different structural "events" (e.g. depleted, emerged, 
split) followed by finding associations between those and SLE 
disease specific traits, like disease activity or severe 
manifestations (e.g. nephritis). 


## Description of the code<a name="code"></a>
The analysis has been performed in R. All the code used to define, 
detect and anaylise DCEs can be found inside the [scripts](https://github.com/vntasis/SLE_spatial_gene_expression/blob/master/scripts/)
folder.
## Cite<a name="cite"></a>
