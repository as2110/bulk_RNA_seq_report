##Deseq2 inputs:
check_kallisto()
##filepath to gene counts tsv/tsv.gz

raw_counts <- "df_mRNA.tsv"

##filepath to sample meta data - minimum 2 columns where column 1 = colheads of gene counts. ncol in gene counts == nrow in sample meta data.
sample_meta_data <- list.files(pattern = "*.tsv", path= ".", full.names = TRUE) %>% stringr::str_subset("design") %>% stringr::str_subset("kallisto")

##column heading for design model for deseq2
deseq2_design_condition <- stringr::str_split(file_path_sans_ext(basename(sample_meta_data)), pattern = "_")[[1]][5]

##contol variable name from the design column
Control <- stringr::str_split(file_path_sans_ext(basename(sample_meta_data)), pattern = "_")[[1]][3]

##Test/model variable name from the design column
Test <- stringr::str_split(file_path_sans_ext(basename(sample_meta_data)), pattern = "_")[[1]][2]

##column heading for column that indicate sample batches if relevant for deseq2
Batch_column_name <- NA

##which stats test to use in deseq2 "Wald" or "LRT"
deseq2_stats_test <- stringr::str_split(file_path_sans_ext(basename(sample_meta_data)), pattern = "_")[[1]][4]

##title for all images/plots/tables.
Plot_title <- stringr::str_split(file_path_sans_ext(basename(sample_meta_data)), pattern = "_")[[1]][6]


##log2 fold change cut off for Volcano plot.
FC_VP <- 2

##pvalue cut off for Volcano plot.
pvalue_VP <- 1e-8

##GSEA inputs
##choose MSigDB pathways: can be "Hallmark", "Kegg" or "Reactome"
pathways <- "Hallmark"


#ORA

##fold change cut of over representation analysis. Values <0 will make list at under expressed genes and >0 will overexpressed genes
FC_ORA <- 1

##p value cut off for over representation analysis.
pvalue_ORA <- 0.05



#### prep raw data

prepare_files(raw_counts, sample_meta_data)


### make deseq2 model

DeSeq2_return <- Deseq2_results(counts = raw_counts, samples = sample_meta_data, design = deseq2_design_condition, control = Control, test = Test, batch = Batch_column_name, stats = deseq2_stats_test)



### annotate results

Annotate_genes_results(res = DeSeq2_return$res,
                       title = Plot_title,
                       control = Control,
                       test = Test,
                       stats = deseq2_stats_test
)
