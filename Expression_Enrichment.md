# This script is for running permutation tests for gene expression enrichment within a subset of genes against the genomic background.


### The function below requires three parameters:
#### 1: query_H0 - the list of genes to query *in this case, this would be all the genes in one inversion*

query_H0 should be a data frame object with a single column, each row representing one gene:

`#head(query_H0)

#                     V1
#1 Ha412HOChr16g00049975
#2 Ha412HOChr16g00047516
#3 Ha412HOChr16g00049976
#4 Ha412HOChr16g00047517
#5 Ha412HOChr16g00047518
#6 Ha412HOChr16g00047519`


#### 2: inv_library - the expression library you are querying against *in this case, average expression of all genes in 10 tissues*

inv_library should be a data frame object with the first column containing gene names with the colname "Gene." All other columns should contain expression values for all genes in the genome (alternatively, if you want to query against a smaller library, include only genes you want to query against):

`#head(inv_library)

#                        Gene Bract Corolla Pollen Ligule Stamen Stem Style        Leaf       Root  Ovary
#1 Ha412HOChr00c00007g00053165   199      27     30    309     27  326    75 113.0000000  60.500000 126.75
#2 Ha412HOChr00c00008g00053166   153       2      5     84     22  532    15  47.7777778   1.333333  24.50
#3 Ha412HOChr00c00014g00053167   553      56     22    666     51  578   210 203.1111111 185.666667 274.00
#4 Ha412HOChr00c00020g00053168     0       0    595     14     61    0     4   0.3333333   0.000000 152.25
#5 Ha412HOChr00c00023g00053169   362     298    113   2163    594  349  1749  32.4444444  66.666667 657.75
#6 Ha412HOChr00c00094g00053170   266     101     10    422     96  530   194  67.6666667 159.166667 161.00`

#### 3: permutations - the number of times to permute the data *anywhere from 10,000 to 100,000 permutations is considered an acceptable balance between accuracy and computational time* 

`permuations <- 10000`

The code follows:

`query_permutes <- function(query_H0, inv_library, permutations){
  #rename query_H0 column
  colnames(query_H0)[1] <- "Gene"
  #extract expression values for query genes
  query_H0 <- merge(query_H0, inv_library, by = "Gene", all.x = T, all.y = F)
  
  #limit the number of genes to search per permutation to 75 to limit computational strain, otherwise use 1/4 of the query genes per permutation
  if(floor(nrow(query_H0)/4) > 75){
    query_temp <- 75
  } else { query_temp <- floor(nrow(query_H0)/4) }
  
  #extract 500 genes per permutation from the library:
  lib_temp <- 500
  
  #if the number of genes extracted from the query is greater than 2, proceed with the analysis. Otherwise, it will not run
  if(query_temp > 2){
    
    colnames_vec <- c()
    
    #run permutations across all columns of the expression library. In our case, this means 10,000 permuations for each tissue
    
    for(column in 2:ncol(inv_library)){
      count = 1 #counter for number of permuations
      enrich_p = 0 #counter for number of times that library < query
      deplete_p = 0 #counter for number of times that library > query
      working_count = 0 #counter for number of times either enrich_p or deplete_p go up (excludes TIES)
      
      #for each permutation:
      while(count <= permutations){
        inv <- dplyr::sample_n(query_H0, size = query_temp, replace = FALSE) #sample query_temp genes from the query library
        ctl <- dplyr::sample_n(inv_library, size = lib_temp, replace = FALSE) #sample 500 genes from the genomic library
        
        #find the average expression value from each:
        inv_mean <- as.numeric(colMeans(inv[column]))
        ctl_mean <- as.numeric(colMeans(ctl[column]))
        
        
        if(ctl_mean > inv_mean){
          deplete_p = deplete_p + 1 # p_deplete increases by 1 if the genomic library average was higher
          working_count = working_count + 1
        }else if(ctl_mean < inv_mean){
          enrich_p = enrich_p + 1 # p_enrich increases by 1 if the query library average was higher
          working_count = working_count + 1
        }
        count = count + 1
      }
      
      #p-value of enrichment is (1-number of enriched permuations)/total number of permuations 
      enrichment_p <- deplete_p/working_count
      depletion_p <- enrich_p/working_count
      
      #output a dataframe with one row per permutation test:
      #colnames_vec = (colnames_vec, colnames(inv)[column])
      temp_df <- data.frame(colnames(inv)[column], enrichment_p, depletion_p)
      rownames(temp_df) <- c()
      colnames(temp_df) <- c("Variable", "p_Enrich", "p_Deplete")
      
      
      if(exists("output_df") == F){
        output_df <- temp_df
      }else{output_df <- rbind(temp_df, output_df)}
    }
    
    return(output_df)
  } else {
    return(NA)
  }
  
}`












