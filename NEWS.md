# SmartPhosExplorer update notes

## Version: 1.0.4

Date: 2024-03-27

CHANGES IN VERSION 1.0.4

------------------------------------------------

Time series clustering:    

* Introduce option for adding 0 timepoint for treatments for clustering and visualization. It simply copies the zero timepoint samples from the selected control. 
* Set seed for making clustering results reproducible.

Preprocessing:

* Along with the multiAssayExperiment object, user can now save the summarized experiment object generated after clicking on the process button in the preprocessing panel. If the user has selected annotations for filtering in the displayed table, then a summarized experiment object with filtered values will be saved. 

-------------------------------------------------------

## Version: 1.0.3

Date: 2024-02-23

CHANGES IN VERSION 1.0.3

------------------------------------------------

Log Info:  

* Log info tab is added to the SmartPhos explorer app to keep track of all the user inputs.
* All the log info can be download as a TSV file.

Differential expression analysis:    

* Volcano plot is made interactive. When the user clicks on the row of differential expression genes table, the corresponding point on volcano plot is highlighted as a star symbol.
* Also, when user click on any colored point on volcano plot, the point turns into star and the corresponding row in the table is highlighted. 

Preprocessing:

* User can subset data using the filtering options present on the top of the table being displayed.
* Option to correct for batch effects is added in the preprocessing step. 

-------------------------------------------------------

## Version: 1.0.2

Date: 2024-01-16

CHANGES IN VERSION 1.0.2

------------------------------------------------

Data preprocessing:  

* SmartPhos app does not require the upload of complete proteomic and phosphoproteomic search of both enriched and un-enriched samples. For examples, users can upload only phospho search on the enriched samples or proteomic search on the unenriched samples for the analyses.   

* A bug in normalization with log2 transformation has been fixed.   

Differential expression analysis:    

* Introduce Volcano plot for the differential expression analysis result.
* User can click on the coloured data points of the volcano plot to get the box plot. This box plot is same as the one generated when clicking on the row of DE table with same ID.

-------------------------------------------------------

## Version: 1.0.1

Date: 2023-12-07

CHANGES IN VERSION 1.0.1

------------------------------------------------

Differential expression analysis: 

* Change the way samples are selected for more flexibility. Users can either choose by sample ID or selecting the applicable treatment(s) and time point(s) for the reference and target groups. The number of samples in each group is updated automatically and a warning is given if a sample appears in both groups.
* Remove dependence on the 'timepoint' column in 'fileTable.txt' to run the analysis.

Time series clustering:

* Uses the mean instead of median when computing a representative statistic across replicates.
* Remove dependence on the replicate column when computing logFC (previously an error is shown if both the 'subjectID' and 'replicate' or 'rep' columns are absent in the 'fileTable.txt').

Enrichment analysis:

* Add phospho-signature enrichment to find enrichment in a site-centric database.
* For the said database: use the PTM signature database version 2.0.0 from https://proteomics.broadapps.org/ptmsigdb/. For more details on the curation of this database, see Krug et al., 2019: https://doi.org/10.1074/mcp.TIR118.000943. In this database, only phosphorylation sites are considered and signature sets whose names start with "KINASE" are removed (since it would be similar to the kinase activity analysis).
* For the algorithm in phospho-signature enrichment: use PTM Signature Enrichment Analysis (PTM-SEA) for Differential Expression and Fisher's exact test for Time series clustering. PTM-SEA was modified from the work of Krug et al.

Kinase activity analysis:

* Add option to adjust the number of permutations to calculate the null distribution.
* Add Fast Gene Set Enrichment Analysis (FGSEA) as an option to estimate the association of kinases to a time series cluster. FGSEA is done by the 'run_fgsea()' function from decoupleR using probabilities that the phospho sites belong to the selected cluster to rank the sites.
* Add option to use FDR

BUG FIXES:

* Fix spline filtering for time series clustering.
* Changes titles of boxplots and time course plots to be the name of the selected phospho site instead of gene.
* Remove dependence on the package 'genefilter'
* Fix color display for the horizonal barplot in kinase activity analysis.
* Fix x-axis display for the heatmap in kinase activity analysis.

-------------------------------------------------------

## Version: 1.0.0

Date: 2023-11-10

------------------------------------------------

First stable release on EMBL server
