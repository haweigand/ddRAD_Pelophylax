# ddRAD_Pelophylax
## Code for ddRAD analysis of Pelophylax water frog data


## Parameter Test
1) Run Stacks with different settings

2) Export vcf-file for each setting

3) Filter loci and individuals fullfilling certain quality criteria 
   
   **filter_vcf.py**
  * input: genepop file + files to filter the file 
  * filter options (e.g. minimal number of individuals/snps, minor allele frequency)
  * output: genepop, fasta, lfmm (lea)

4) Generate summary statistics from datasets with different qualtiy criteria
   
   **summary_statistics.py**
  * script to generate summary statistic for parameter test based on replicated individuals 
  * input: folder with genepop files (output from filter_vcf.py)
  * settings: thresh => summary of the dataset: proportion of individuals with more SNP-calling errores than the threshold in a certain setting
  * output: 
     * statistic_genepop_individuals.txt => summary statistics for each individual
     * statistic_genepop_summary.txt => summary statstics for the total dataset

## Final analysis
1) Run stacks with parameters chosen by the parameter test & export vcf file

2) Filter loci and inidviduals with the filtering settings choosen from the parameter test
   **filter_vcf.py**
   (see above)

3) Run a PCA in Lea

4) Use PCA-Result to identify individuals for diagnostic loci identification
  **selectIndFromPCA.py**
  * input
     * lea pca projections
     * genepop file
     * list with names of individuals used in pca (lfmm_ind)
     * boarders of pca clusters
  * setting: propOfInd => use individuals that are centrally located in the PCA clusters
  * output: list of individuals close to center of pca cluster with lowest proportion of missing data

5) Filter vcf file for diagnostic loci
  **getDiagnosticLoci.py**


 
   
