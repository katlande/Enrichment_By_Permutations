A function that queries if a set of genes have the same TAIR homologs more often than expected given the genomic frequency of genes with that TAIR homolog

````r

#### Input files ###

#annot = annotation file with Homolog (TAIR homolog determined through blast), Gene (Ha412HOv2.1 gene annotation), and Annotation (TAIR homolog annotation from TAIR) columns.

#Homolog                  Gene        Annotation
#1 AT1G01020 Ha412HOChr14g00040714    ARV1 family protein;(source:Araport11)
#2 AT1G01060 Ha412HOChr15g00043361    LHY encodes a myb-related putative transcript


#Extract all TAIR genes that are homologous to at least 3 Ha412HOv2.1 genes:
t<-as.data.frame(table(annots$Homolog))
t<-as.character(t$Var1[t$Freq>=3])
t<-annots[annots$Homolog%in%t,]


#Given all genes inside the inversion of interest as a character vector, inv_genes, extract all TAIR homologs annotated to genes in the inversion that are also annotated to at least 3 genes in the Ha412HOv2.1 reference genome:

homologs2query <- annots$Homolog[annots$Gene%in%inv_genes & annots$Homolog%in%t$Homolog] #character vector

#Function to find inv and non-inv frequencies of each homolog in a set of inverted genes:

#inputs:
#homologs2query (see above)
#inv_genes (see above)
#statistic = "CHI" or "FISHER" depending on which test you want to run.
trawl_homologs <- function(homologs2query, inv_genes, statistic="CHI"){
  
  
  #columns for the output:
  Homolog <- c()
  Inverted <- c()
  Non_Inverted <- c()
  Expected <- c()
  p <- c()
  
  inv_genesl <- length(inv_genes)
  non_inv_genes <- nrow(gff)-inv_genesl
  
  #loop to extract data for each homolog:
  for(homolog in homologs2query){
    
    #find all the genes homologous to each TAIR homolog:
    h_data <- annots[annots$Homolog == homolog,]
    inv_homolog <- length(h_data$Gene[h_data$Gene%in%inv_genes])
    non_inv_homolog <- nrow(h_data) -inv_homolog
    
    #Contingency table::
    #_______________#______________#___________________#
    #               #    homolog   #    other genes    #
    #_______________#______________#___________________#
    #   Inverted    #1       ?     #3        ?         #
    #_______________#______________#___________________#
    # non-Inverted  #2      ?      #4        ?         #
    #_______________#______________#___________________#
    #
    #1: inv_homolog
    #2: non_inv_homolog
    #3:
    other_inv <- (inv_genesl - inv_homolog)
    #4:
    other_non_inv <- (non_inv_genes - non_inv_homolog)
    
    #actual inv: inv_homolog
    #actual non-inv: non_inv_homolog
    #expected inv:
    exp_inv <- (inv_genesl/nrow(gff))*nrow(h_data)
    #expected non-inv:
    exp_non_inv <- (non_inv_genes/nrow(gff))*nrow(h_data)
    
    #if more homologs are expected to be inside the inversion, skip this gene
    if(exp_inv >= inv_homolog){
      next
    }
    
    #define statistical test:
    if(statistic=="CHI"){
      res <- stats::chisq.test(matrix(c(inv_homolog, non_inv_homolog, other_inv, other_non_inv),nrow=2))
      chi_p <- res$p.value
      p <- c(p, chi_p)
    } else if(statistic=="FISHER"){
      res <- stats::fisher.test(matrix(c(inv_homolog, non_inv_homolog, other_inv, other_non_inv),nrow=2), y = NULL, alternative = "greater", conf.level = 0.95)
      fisher_p <- res$p.value
      p <- c(p, fisher_p)
    }
    
    Homolog <- c(Homolog, homolog)
    Inverted <- c(Inverted, inv_homolog)
    Non_Inverted <- c(Non_Inverted, non_inv_homolog)
    Expected <- c(Expected, exp_inv)
    
  }
  
  if(statistic=="CHI"){
    output <- data.frame("Homolog"=Homolog, "Inverted"=Inverted, "Non_Inverted"=Non_Inverted, "Expected"=Expected, "ChiSq_p"=p)
  }else if(statistic=="FISHER"){
    output <- data.frame("Homolog"=Homolog, "Inverted"=Inverted, "Non_Inverted"=Non_Inverted, "Expected"=Expected, "Fisher_p"=p)
  }
  
  rownames(output) <- c()
  output <- output[output$Inverted > 1,]
  
  return(output)
}

#Function outputs enrichment values for each TAIR homolog (enriched inside inversion compared to the full genome), in a df with the following column names:
#Homolog
#Inverted (# copies inside inversion)
#Non_Inverted (# copies outside inversion)
#Expected (expected # copies inside the inversion)
#p (raw p-value)

trawled_homologs <- trawl_homologs(homologs2query, inv_genes, statistic="CHI")
trawled_homologs$FDR <- p.adjust(trawled_homologs$p, method = "fdr")

#Querying the background rate of gene clusters with the same TAIR homolog in non-inverted sections of the genome:

#Ratio of enriched to non-enriched homologs:
inv_Ratio <- ((nrow(Final_Output[Final_Output$FDR <= 0.01,]) - nrow(Final_Output[Final_Output$FDR > 0.01,])) / nrow(Final_Output))


#Create non-inverted regions to trawl:

chr_st <- c()
chr_e <- c()
Chromosomes <- c("Ha412HOChr02", "Ha412HOChr04", "Ha412HOChr08")#Chromosomes without inversions
for(chr in Chromosomes){
  min <- min(gff$V4[gff$V1 == chr])
  max <- max(gff$V5[gff$V1 == chr])
  chr_st <- c(chr_st, min)
  chr_e <- c(chr_e, max)
}
chr_data <- data.frame("Chr"=Chromosomes, "Start"=chr_st, "End"=chr_e)

#chr_data is a df with the bounds of three inversion-less chromosomes

Query_Size = 10000000 #number of bp to query as a "pseudo-inversion"
permutations = 0
higher = 0
lower = 0
neutral = 0

#query non-inverted "pseudo-inversions" 10,000 times and count how often they are more enriched for TAIR homologs than the average for inverted regions:

while(permutations < 10000){
  
  #randomly choose a region of the non-inverted windows to query:
  temp_row <- sample(1:3, 1, replace=F)#chrs without invs only
  min <- chr_data$Start[temp_row]
  max <- chr_data$End[temp_row]-Query_Size
  temp_window <- sample(min:max, 1, replace=F)
  start <- temp_window
  end <- temp_window + Query_Size
  
  #extract all the genes in that window from the gff (df with no colnames):
  temp_genes <- gff$V9[gff$V1 == as.character(chr_data$Chr[temp_row]) & gff$V4 <= end & gff$V5 >= start]
  homologs2query <- t[t$Gene%in%temp_genes,]
  homologs2query <- unique(homologs2query$Homolog)
  
  #trawl homologs in the pseudo-inversion:
  output <- trawl_homologs(homologs2query, temp_genes, statistic="CHI")
  output$FDR <- p.adjust(output$ChiSq_p, method = "fdr")
  
  #if there are enriched homologs in the query window:
  if(nrow(output) > 1){
    ratio <- ((nrow(output[output$FDR <= 0.01,]) - nrow(output[output$FDR > 0.01,])) / nrow(output))
    
    #check if the frequency of clusters is higher in the random non-inverted query or the whole genome:
    if(ratio > inv_Ratio){higher = higher + 1} else if(ratio < inv_Ratio){lower=lower+1} else { neutral = neutral + 1}
  } else { neutral = neutral + 1}#if the values are equal
  
  permutations = permutations + 1
}

enrichment_p <- lower/(permutations-neutral)
depletion_p <- higher/(permutations-neutral)


````