tSNE_1:
Simple workflow to produce tSNEs using data from flow cytometry

Data Origin:
Origin of data to be used in this analysis: https://flowrepository.org/public_experiment_representations/2560

Only data from sample 1 was used. sample 1 day 8 & and day 11 (below)

1) CYTOF_DS_SPleen_day8_PBS_1.fcs
2) CYTOF_DS_Spleen_day11_PBS_1.fcs

Prepare files for R:
1. Convert .fcs files to .csv files. 

    a) This can be done in R or using other software. I used the website floreada.io

    b) Concatenate .csv files. The files should be concatenated before doing a tSNE so that 
        a single tSNE can be generated.

2. If analyzing a specific phenotype within the tSNE, the original data does not tag 
   specific phenotypes and this must be done in R or separate software. The file "sample_subset.csv"
   is the subset from the original .fcs file that was created for this demonstration.

    a) If desired, create a new file within floreada.io (e.g. select populations that are CD45(+) 
        CD3(+) CD4(+)), and the events can be exported and concatenated with non-duplicate events. 
        Use BASH to assign a unique ID to ID later within R.

3. Import data into R and run the tSNE
4. Follow comments within the provided Rscript

