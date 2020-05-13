# Content

- DCEs: Contains 4 bed files, one per group studied, with the coordinates of the identified DCEs. Each bed file has 4 columns.
There is one column reporting the chromosome followed by two columns associated with the start and end coordinates repspectively.
The fourth column reports the average bin-signal within the respective domain.
- DEGs: Includes bed files reportng the statistically significant differentially expressed genes (DEGs) of each patient group
respectively to the healthy group, as identified by our analysis. There are three files per patient group,
one that reports every DEG (total), one that reports only the overexpressed genes (up), and one that reports only the
underexpressed genes (down). All of them list DEG coordinates. Moreover, the fourth column reports the
log<sub>2</sub>(Fold Change) value and the fifth column indicates the name of the gene.
- DCE_classes: Involves files that describe the spatial re-organization of DCEs,
observed when we compared patient gourps with the healthy group. The files are 
structured in the same way as described above for the "DCEs" folder. Furthermore, in some files both versions of a DCE,
that one identified in the healthy group (H) and the corresponding one in the patient group (D) are reported, followed by the 
Jaccard index value that indicates the overlap between those two.
- Disruptors: Catalogues of Disruptor genes per patient group as identified by our analysis.
- WGCNA_modules: Bed files reporting the members of the WGCNA modules as identified by our analysis. Gene coordinates, strand
and gene name are listed.
