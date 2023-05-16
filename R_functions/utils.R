#BiocManager::install(c("multiGSEA", "fgsea", "clusterProfiler",
#                       "enrichplot", "GOSemSim", "goseq"), ask=FALSE)
#BiocManager::install(c("PCAtools"), ask = FALSE)
###load libraries ####

# if (!require("librarian")) {install.packages("librarian")}
# librarian::shelf(magrittr, tidyverse, ggpmisc, ggpubr, extrafont, fs, tools, ggplotify, grid, kableExtra, RColorBrewer)
# librarian::shelf(yaml, rhdf5, tximport, DESeq2, apeglm)
# librarian::shelf(AnnotationDbi, org.Hs.eg.db, biomaRt, EnsDb.Hsapiens.v86, TxDb.Hsapiens.UCSC.hg38.knownGene)
# librarian::shelf(EnhancedVolcano, PCAtools, pheatmap, ComplexHeatmap)
# librarian::shelf(msigdbr, fgsea, clusterProfiler, goseq, GOSemSim, enrichplot, enrichR)
# librarian::shelf(Pi, "hfang-bristol/XGR")
# librarian::shelf(WGCNA, CEMiTool, GWENA, BioNERO, corto, KBoost, "jpvert/tigress", lionessR, RTN)
# librarian::shelf(httr, jsonlite, dorothea, decoupleR, TFEA.ChIP, CeTF, RcisTarget, RcisTarget.hg19.motifDBs.cisbpOnly.500bp)

#multiGSEA
##general plots
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(extrafont)
library(fs)
library(tools)
library(ggplotify)
library(grid)
library(kableExtra)
library(RColorBrewer)

##differential expression

library(yaml)
library(tximport)
library(DESeq2)
library(apeglm)

#nice plots
library(EnhancedVolcano)
library(PCAtools)
library(pheatmap)
library(ComplexHeatmap)

#Annotations
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library("EnsDb.Hsapiens.v86")

#BiocManager::install(c("EnsDb.Hsapiens.v86", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
##GSEA
library(msigdbr)
library(fgsea)
#library(multiGSEA)
#library(EGSEA)

#ORA
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(enrichplot)
library(goseq)
library(GOSemSim)
library(enrichR)

## other GSEA
library(Pi)
library(XGR)

##GRN
library(WGCNA)
library(CEMiTool)
library(BioNERO)
library(GWENA)
library(corto)

#TF analysis
library(httr)
library(jsonlite)
library("dorothea")
library("decoupleR")
library(TFEA.ChIP)
library(CeTF)
library(RcisTarget)
library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)

# setClass("DeSeq2_return",
#          representation(dds= "DESeqDataSet",
#                         de_seq= "DESeqDataSet",
#                         res ="DESeqResults"
#                         ))

### kallisto ####

check_kallisto <- function(){
params <- read_yaml("config.yml")

if(params$kallisto){
  print("kallisto input detected")
  #yml_site_opts(output_dir = file.path("RNA_seq_reports/kallisto"))

  #design_files <- list.files(path='designs_kallisto/', pattern = "design_")
  design_dir <- "./design_files/design_kallisto"

  # load in the sample metadata file
  sample_meta_data <<- list.files(pattern = "*.tsv", path = design_dir, full.names = TRUE) %>% stringr::str_subset("design") %>% stringr::str_subset("kallisto")
  sample_table <- read.delim(sample_meta_data) %>% dplyr::rename("Sample" = 1)


  if (dir.exists(params$kallisto_dir)) {
    path <- file.path(params$kallisto_dir, sample_table$Sample, "abundance.h5")
    sample_table <- dplyr::mutate(sample_table, path = path)

    if(params$species == "human"){
      # Create the tx2gene file that will map transcripts to genes
      Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id", "symbol"))
      Tx <- as_tibble(Tx)
      Tx <- dplyr::rename(Tx, target_id = tx_id, ens_gene = gene_id, ext_gene = symbol)
      Tx <- dplyr::select(Tx, "target_id", "ens_gene", "ext_gene")

    }else if(params$species == "mouse"){
      Tx <- transcripts(EnsDb.Mmusculus.v79, columns=c("tx_id", "gene_id", "symbol"))
      Tx <- as_tibble(Tx)
      Tx <- dplyr::rename(Tx, target_id = tx_id, ens_gene = gene_id, ext_gene = symbol)
      Tx <- dplyr::select(Tx, "target_id", "ens_gene", "ext_gene")
    }else if(params$species == "macaque"){
      orgSymbols <- keys(org.Mmu.eg.db, keytype="ENSEMBL")
      mart <- useMart(dataset = "mmulatta_gene_ensembl", biomart='ensembl')
      tx2gene <- getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version','entrezgene_id'),
                       filters = 'ensembl_gene_id',
                       values = orgSymbols,
                       mart = mart)
      Tx <- tx2gene %>%   dplyr::select(ensembl_transcript_id, ensembl_gene_id)
      colnames(Tx) <- c("TXNAME", "GENEID")

    }else if(params$species == "pig"){
      ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="sscrofa_gene_ensembl",
                        host="uswest.ensembl.org",
                        ensemblRedirect = FALSE)
      orgSymbols <- unlist(getBM(attributes = 'ensembl_gene_id', mart=ensembl))
      mart <- useMart(dataset = "sscrofa_gene_ensembl", biomart='ensembl', host="uswest.ensembl.org")
      tx2gene <- getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version','entrezgene_id'),
                       filters = 'ensembl_gene_id',
                       values = orgSymbols,
                       mart = mart)
      Tx <- tx2gene %>%   dplyr::select(ensembl_transcript_id, ensembl_gene_id)
      colnames(Tx) <- c("TXNAME", "GENEID")
    }else if(params$species == "rabbit"){
      ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="ocuniculus_gene_ensembl",
                        host="uswest.ensembl.org",
                        ensemblRedirect = FALSE)
      orgSymbols <- unlist(getBM(attributes = 'ensembl_gene_id', mart=ensembl))
      mart <- useMart(dataset = "ocuniculus_gene_ensembl", biomart='ensembl', host="uswest.ensembl.org")
      tx2gene <- getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id', 'ensembl_transcript_id_version','entrezgene_id'),
                       filters = 'ensembl_gene_id',
                       values = orgSymbols,
                       mart = mart)
      Tx <- tx2gene %>%   dplyr::select(ensembl_transcript_id, ensembl_gene_id)
      colnames(Tx) <- c("TXNAME", "GENEID")
    }else{ print('please supply a valid species within the config.yml file')}

    # import Kallisto transcript counts into R using Tximport
    Txi_gene <-  tximport(path,
                          type = "kallisto",
                          tx2gene = Tx,
                          txOut = FALSE, # TRUE outputs transcripts, FALSE outputs gene-level data
                          countsFromAbundance = "lengthScaledTPM",
                          ignoreTxVersion = TRUE)
    # Write the counts to an object (used for metadata and clustering)
    df_mRNA <- Txi_gene$counts %>%
      round() %>%
      data.frame()
    colnames(df_mRNA) <- sample_table$Sample

    meta_data <- sample_table
    rownames(meta_data) <- meta_data$Sample

    # if (file.exists(paste0('designs_kallisto/',design_files[1]))) {
    #   for (i in design_files){
    #   meta_data <- read.table(paste0('designs_kallisto/',i), sep=",", header = TRUE)
    #   rownames(meta_data) <- meta_data$Sample
    #   df_mRNA_tmp = df_mRNA[,rownames(meta_data)]
    #   all(rownames(meta_data) %in% colnames(df_mRNA_tmp))
    #   assign(paste("meta_data", i, sep = "."), meta_data)}}
    #   else {
    #     print("No design files were detected please add a file called design_<test>_<control>_<test>_<column>.csv. Please refer to documentation on github for more ifnormation")
    #     }

    write.table(df_mRNA %>% rownames_to_column("id"),
                paste0("genes.tsv"),
                col.name=TRUE,
                sep="\t",
                na = "NA",
                row.names=FALSE,
                quote=FALSE)
    }
} else {
  design_dir <- "design_files/design_featurecounts"
  sample_meta_data <<-
    list.files(pattern = "*.tsv",
               path = design_dir,
               full.names = TRUE) %>% stringr::str_subset("design") %>% stringr::str_subset("featurecounts")
}
}
## prep files #####




prepare_files <- function(raw_counts, sample_meta_data){

  ifelse(!dir.exists("results"), dir.create("results"), FALSE)


  ifelse(!dir.exists("raw_data"), dir.create("raw_data"), FALSE)
  ifelse(!dir.exists("meta_data"), dir.create("meta_data"), FALSE)
  file.copy(from = raw_counts, to = "raw_data", overwrite = TRUE)
  file.copy(from = sample_meta_data, to = "meta_data", TRUE)

  raw_counts <<- paste0("raw_data/",path_file(raw_counts))
  sample_meta_data <<- paste0("meta_data/",path_file(sample_meta_data))




}



#### deseq2 ####
Deseq2_results <- function(counts, samples, design, control, test, stats = "Wald", batch = NA){

  #create results directory for this control
  ifelse(!dir.exists(paste0("results/", Plot_title, "_", control, "_vs_", test)), dir.create(paste0("results/", Plot_title, "_", control, "_vs_", test)), FALSE)
  ifelse(!dir.exists(paste0("results/", Plot_title, "_", control, "_vs_", test, "/plots")), dir.create(paste0("results/", Plot_title, "_", control, "_vs_", test, "/plots")), FALSE)
  plots_dir <<- paste0("results/", Plot_title, "_", control, "_vs_", test, "/plots/")
  results_dir <<- paste0("results/", Plot_title, "_", control, "_vs_", test, "/")

  ## read the raw data file
  raw <- read_tsv(counts) %>%
    mutate(
      hgnc_symbol = mapIds(
        EnsDb.Hsapiens.v86,
        keys = id,
        keytype = "GENEID",
        column = "SYMBOL"
      ),
      .before = 1
    ) %>%
    drop_na() %>%
    dplyr::select(-id) %>%
    group_by(hgnc_symbol) %>%
    summarise(across(everything(), ~ sum(.))) %>%
    mutate(
      id = mapIds(
        EnsDb.Hsapiens.v86,
        keys = hgnc_symbol,
        keytype = "SYMBOL",
        column = "GENEID"
      ),
      .before = 1
    ) %>% dplyr::select(-2)

  #View(raw)
  ## transcript is a list containing all the transcript names

  transcripts <- pull(raw, 1)


  ## Create a matrix from the tsv file. cts is a matrix with numerical data only for all 30 samples

  cts <- as.matrix(raw[,-1])

  ##annotate row names in cts with the pulled genes names vector.

  rownames(cts) <- transcripts

  ##read in the meta data file
  sampleinfo <- read.delim(samples) %>% rowid_to_column()

  sampleinfo <- sampleinfo[, c(2,1,3:ncol(sampleinfo))]



  des <- paste0("~", design)

  if (!is.na(batch)) {
    sampleinfo <- mutate(sampleinfo, batch = as.factor(sampleinfo[[batch]])
                         )
    des <- paste0("~", batch, "+", design)
  }




  if (nrow(sampleinfo) == ncol(cts)) {
     dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = sampleinfo,
                                design = as.formula(des))

  ##Set control to reference level
  dds[[design]] <- relevel(dds[[design]], ref = control)
  dds[[1]] <- as.factor(dds[[1]])

  }

  ## filter for expressed genes.

  keep <- rowSums(assay(dds) >= 10) >= round((nrow(sampleinfo)/5))
  #keep <- is_expressed >= 5
  #table(keep)
  dds <- dds[keep,]


  if (stats == "LRT") {
    de_seq <- DESeq(dds, test = stats, reduced = ifelse(is.na(batch),
                                                         as.formula("1"),
                                                         as.formula(batch)
                                                         ))



  }

  if (stats == "Wald") {
     de_seq <- DESeq(dds)
  }


  ##Change the name= from resultNames(de_seq) depending on the chosen design.

  contrast_list <- DESeq2::resultsNames(de_seq)
  contrasts <- contrast_list[str_which(contrast_list, test)]

  ##Run results
  res <- DESeq2::results(de_seq,
                         name = contrasts,
                         tidy = FALSE,
                         independentFiltering = TRUE)



  #return(new("DeSeq2_return",
   #           dds= dds,
    #          de_seq = de_seq,
     #         res = res
      #        ))


  vsd <- DESeq2::vst(dds,blind=TRUE)

  DeSeq2_return <-list("dds" = dds,
                       "vsd" = vsd,
                       "de_seq" = de_seq,
                       "res" = res,
                       "contrast_list" = contrast_list)

  return(DeSeq2_return)


}

##Visualising library size ####
Library_size_check <- function(title, dds, fill = NA, design){

  cts <- assay(dds)
  total_reads <- colSums(cts)/1000000
  sampleinfo <- colData(dds) %>% as.data.frame() %>%
                dplyr::mutate(Millions_of_Reads = total_reads)
  cols <- colnames(sampleinfo)







  Bar <- ggbarplot(sampleinfo, x = "rowid", y = "Millions_of_Reads",
                   fill = ifelse(is.na(fill), design, fill), ## change fill color by  Batch or other column
                   palette = c("dark red", "dark green", "navy", "black"),
                   color = "black",            # Set bar border colors to white
                   # palette = "jco",            # jco journal color palette. see ?ggpar
                   # sort.val = "desc",          # Sort the value in descending order
                   sort.by.groups = FALSE,     # Don't sort inside each group
                   label = sampleinfo[[design]],
                   lab.size = 1.5,
                   lab.pos = c("out"),
                   title = title,
                   submain = "Library Size",
                   xticks.by = 1,
                   legend = "bottom",
                   legend.title = ifelse(is.na(fill), design, fill),
                   xlab = "Sample",
                   ylab = "Millions of Reads"


  ) %>% ggpar(ggtheme = theme_pubr(),
              font.title = c("bold", "brown", 18),
              font.subtitle = c("bold", "dark grey", 14),
              font.caption = c("bold.italic", "royal blue", 12),
              font.legend = c("bold", "black", 8),
              font.x = c("bold", "black", "11"),
              font.y = c("bold", "black", "11"),
              font.tickslab = c("bold", "black", "8"),
              font.family = "Calibri"
  )





  #Visualize bar chart
  return(Bar)

  ##to save use:
  #ggsave(paste0(plots_dir, title, "_library_size_BarPlot.png"))

}



Library_VP <- function(vsd, title, design, fill = NA, test, control) {

  ##convert matrix into table for VP plots
  Log2_library_data <- assay(vsd) %>% as.table() %>% as.data.frame()
  sampleinfo <- as.data.frame(colData(vsd))
  colnames(Log2_library_data) <- c("TranscriptID", colnames(sampleinfo)[1], "Log2_Counts_per_million")


  Log2_library_data <- Log2_library_data %>%
    inner_join(sampleinfo, by=colnames(sampleinfo)[1])

  ###Box and Violin plots using ggpubr format.

  ##To select specfic samples to plot in violin plot
  pattern <- c(control, test)
  fill <- ifelse(is.na(fill), design, fill)
  violin_selection <- pull(dplyr::filter(Log2_library_data, Log2_library_data[ , design] %in% pattern), rowid)

  ##Plot Violins
  boxp <- (ggviolin(Log2_library_data,
                    x = "rowid",
                    #x = "SampleName",
                    y = "Log2_Counts_per_million",
                    add = c("median_iqr", "boxplot"),
                    add.params = list(fill = "white", color = "black"),
                    draw_quantiles = 0.5,
                    palette = c("dark red", "dark blue", "dark green", "black", "orange"),
                    fill = fill,   ##chose what factor to color by
                    select = violin_selection,
                    facet.by = design,
                    #label = colnames(sampleinfo)[1],


          )) %>%


    ggpar(ylab = bquote(~"Log"[2]~"Counts per million"),
          xlab = "Sample",
          title = title,
          #subtitle = contrasts,
          legend.title = fill,
          #caption = "Normalised Library distribution",
          legend = "top",
          #legend.title = list(fill = "Batch"),
          font.title = c("bold", "brown", 18),
          font.subtitle = c("bold", "dark grey", 14),
          font.caption = c("bold.italic", "royal blue", 12),
          font.legend = c("bold", "black", 8),
          font.x = c("bold", "black", "11"),
          font.y = c("bold", "black", "11"),
          font.tickslab = c("bold", "black", "8"),
          font.family = "Calibri",
          x.text.angle = 0,


          ggtheme = theme_pubr()  ##can use other inbuilt themes or customise

    )

  return(boxp)
  #ggsave(paste0(plots_dir, title, "_library_size_plots.png"))
}


## annotate ####
Annotate_genes_results <- function(res, title, control, test, stats){

##Use AnnotationDbi and org.Hs.eg.db to get ENTREZIDs and biomart() gene symbols because this gives the least missing values.
## to get names of keys
# keytypes(org.Hs.eg.db)
## to see what you can retrieve
# columns(org.Hs.eg.db)
##change org database as needed.
##res is the results from results(de_seq).


  res$GENENAME <- mapIds(org.Hs.eg.db,
                         #EnsDb.Hsapiens.v86,
                         keys=row.names(res),
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")

  res$ENTREZID <- mapIds(#org.Hs.eg.db,
                         EnsDb.Hsapiens.v86,
                         keys=row.names(res),
                         column="ENTREZID",
                         keytype="GENEID",
                         multiVals="first")

  res$hgnc_symbol <- mapIds(#org.Hs.eg.db,
    EnsDb.Hsapiens.v86,
    keys=row.names(res),
    column="SYMBOL",
    keytype="GENEID",
    multiVals="first")


  ###Biomart to get Gene Symbols
  #ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  #mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = ensembl)

  ##convert res to data.frame and merge symbols.
  ##Write to CSV to save the output for later use.

  results_annotated <- as.data.frame(res) %>%
    tibble::rownames_to_column("ensembl_gene_id") #%>%
    #dplyr::left_join(mapping, by = "ensembl_gene_id") %>%
    #dplyr::filter(duplicated(ensembl_gene_id) == FALSE) #%>%
    #na_if("")


  write.table(results_annotated,
              paste0(results_dir, "deseq2_results_annotated.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              row.names=FALSE,
              quote=FALSE)




  ##Viewing a selection of top hits from the Results table using ggtextable()

  #res_selection <- dplyr::filter(results_annotated, padj < 0.05)


  res_selection <- dplyr::select(results_annotated,
                                 ensembl_gene_id,
                                 Gene_Symbol = hgnc_symbol,
                                 Mean_Normalised_Count = baseMean,
                                 test_statistic = stat,
                                 Log2FoldChange = log2FoldChange,
                                 Adjusted_p_value = padj) %>%
    #arrange(desc(abs(test_statistic))) %>%
    arrange(Adjusted_p_value)




   #res_selection <- dplyr::filter(res_selection, Adjusted_p_value < 0.05)



  res_selection <- res_selection[1:30,]  ###Change this line increase rows to display

  write.table(res_selection,
              paste0(results_dir,"top_genes.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              row.names=FALSE,
              quote=FALSE)

  ##Create table with ggtextable

  RESULTS_TABLE <- ggtexttable(res_selection, rows = NULL,
                               cols = c("Ensembl ID", "Gene Symbol", "Mean Normalised Count", paste0(stats, " statistic"), "Log2 Fold Change", "Adjusted P-value"), theme = ttheme("classic", colnames.style = colnames_style(color = "black", face = "bold", size = 10, linecolor = "black", fill = "light grey", linewidth = 1)), ) %>%

    table_cell_font(row = 2:(nrow(res_selection)+1), column = 1:ncol(res_selection), face = "bold", size = 8)


  ##To conditionally color or highlight cells
  # Specify colors for significant log-fold-changes
  res_selection <- res_selection %>%
    dplyr::mutate(
      fill = ifelse(Log2FoldChange > 0, "green", "red"),
      color = "black"
    )

  # Coloring the table conditionally using `ggpubr::table_cell_bg()`

  for(j in 1:nrow(res_selection)){
    row = j+1
    column = which(colnames(res_selection) == "ensembl_gene_id")
    RESULTS_TABLE <- table_cell_bg(
      RESULTS_TABLE, row = row, column = column,
      fill = res_selection$fill[j], color = res_selection$color[j]
    )
  }

  ##Reset to original selection table
  res_selection <- res_selection %>% dplyr::select(-fill, -color)


  # Add title and footnotes
  RESULTS_TABLE <- RESULTS_TABLE %>%
    tab_add_title(text = paste0("Table 1 - Top ", nrow(res_selection), " differentially expressed genes"), face = "bold", size = 12, family = "Calibri") %>%
    tab_add_title(text = paste0(title, " - ", control, " vs ", test) , face = "bold", size = 16, color = "brown", family = "Calibri", padding = unit(0.2, "line"))

  #%>% tab_add_footnote(text = contrasts[i], size = 10, face = "italic", color = "royal blue", )




  ##create ranked file GSEA
  #order values by Wald stat - can also use sign(log2FoldChange) * -log10(pvalue)
  #rnk file has two unnamed columns - gene IDs and ranks. ENTREZID is being used here as most online tools require it.

  gseaInput <-  results_annotated %>%
    dplyr::filter(!is.na(hgnc_symbol), !is.na(stat)) %>% dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
    arrange(stat)


  gsea_input_rnk <- dplyr::select(gseaInput, hgnc_symbol, stat)
  ##return(gsea_input_rnk)
  write.table(gsea_input_rnk, paste0(results_dir, "gsea_input.rnk"),col.name=FALSE,sep="\t",row.names=FALSE,quote=FALSE)

  ##save table image
  RESULTS_TABLE
  ggsave(paste0(plots_dir, title, "_", control, "_vs_", test,  "_top_hits_tables.png"))

  return(RESULTS_TABLE)

}

#DispPLot ####


Dispersion_Plot <- function(de_seq, title, ymin = NA) {

  ###Dispersion plot can be done using in-built function in De-Seq2. It uses Base R to plot.

  #Disperion_plot <- plotDispEsts(de_seq, returnData = FALSE, main = "In built DispEsts Plot function")

  ###To manually plot dispersions using ggplot use code below

  #return(class(de_seq))

  ##Variables

  CV = FALSE     ##set this to TRUE wanted to plot variance instead of dispersion
  f <- if (CV) sqrt else I ##this is allow plotting variance if needed

  ##Extract metadata from de_seq object and set these to x and y axis

  px = mcols(de_seq)$baseMean
  sel = (px>0)
  px = px[sel]

  py = f(mcols(de_seq)$dispGeneEst[sel])

  ifelse(!is.na(ymin), ymin,
         ymin <- 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1))  # ##set y-axis limit


  disp_plot <- (data.frame(baseMean=px, dispGeneEsts=py, name=rownames(de_seq),
                           dispFit = mcols(de_seq)$dispFit,
                           dispOutlier = mcols(de_seq)$dispOutlier,
                           final = f(dispersions(de_seq))) )

  ###In case of setting a lower y-min, the below codes adds a new column
  disp_plot <-  (mutate(disp_plot, dispGeneEstsPlot = ifelse(disp_plot$dispGeneEsts > ymin, disp_plot$dispGeneEsts, ymin)))  %>%


    ###this make the different dispersion counts into a single column.
    pivot_longer(., cols = c(dispGeneEstsPlot,dispFit, final),
                 names_to = "DispersionType",
                 values_to = "DispersionCount") %>%

    ## this is to mark out the outliers
    mutate(shape = ifelse(
      dispGeneEsts <= ymin,
      "Low",
      ifelse(dispOutlier, "Outlier", "Normal")
    ))



  Disp_Plot <- (ggplot(disp_plot,
                       aes(x=baseMean,
                           y=DispersionCount)) +
                  geom_point(aes(color = DispersionType,
                                 shape = shape))
                +scale_color_manual(
                  name = 'Dispersion Estimates',
                  labels = c("dispFit" =  "Fitted", "dispGeneEstsPlot" = "Gene Estimates", "final" = "Final"),
                  values = c("red", "black", "royal blue"))
                +scale_shape_manual(
                  #name = 'Dispersion Estimates',
                  values = c("Low" =  6, "Outlier" =  10, "Normal" = 20))
                +guides(colour = guide_legend("Dispersion Estimates"),
                        shape = "none")
  )  %>%

    ggpar(xscale ="log10",
          yscale = "log10",
          format.scale = TRUE,
          ylim = c(ymin, max(disp_plot$dispGeneEsts)),
          legend.title = list(color = "Dispersion Estimates", shape = "Outliers"),
          title = title,
          #submain = contrasts,
          #caption = "Dispersion Plot",
          ggtheme = theme_pubr(),
          xlab = "Normalised Mean Count",
          ylab = "Dispersion",
          font.title = c("bold", "brown", 18),
          font.subtitle = c("bold", "dark grey", 14),
          font.caption = c("bold.italic", "royal blue", 12),
          font.legend = c("bold", "black", 8),
          font.x = c("bold", "black", "11"),
          font.y = c("bold", "black", "11"),
          font.tickslab = c("bold", "black", "8"),
          font.family = "Comic Sans MS",
    )



  return(Disp_Plot)
  #ggsave(paste0(plots_dir, Dispersion_Plot", title, ".png"))
  ###Need to work out a better way to show outliers and a manual legend for it.




}
#MAPlot ####

MA_Plot <- function(de_seq, test,
                    res,
                    FDR = 0.05,
                    FC = 1,
                    gene_annotations = paste0(results_dir,"deseq2_results_annotated.tsv"),
                    title) {

  ###MA Plots
  ##LFC Shrinkage using apelgm

  ##can easily edit this function so it can use alternate contrasts
  ##just re-run the re DESeq2::results() function with alternate contrast

  contrast_list <- DESeq2::resultsNames(de_seq)
  contrasts <- contrast_list[str_which(contrast_list, test)]
  #res <- DESeq2::results(de_seq, name = contrasts, tidy = FALSE, independentFiltering = TRUE)

  ##read in annotated results table
  results_annotated <- read.delim(gene_annotations) %>%
    as.data.frame()


  resLFC <- lfcShrink(de_seq, coef= contrasts, type="apeglm", res = res)

  ###MA Plot using in-buit function in DE-seq2
  #plotMA <- plotMA(resLFC, ylim=c(-2,2))
  #plotMA

  ###How plot the out of bound data points when ylim < data.


  ###MA Plot manually using ggmaplot

  maplot_apeglm <- ggmaplot(resLFC,
                            fdr = FDR,
                            fc = FC,
                            size = 1,
                            palette = c("#B31B21", "#1465AC", "darkgray"),
                            genenames = results_annotated$hgnc_symbol,
                            #genenames = rownames(resLFC),
                            legend="top", top = 5,
                            font.family = "arial",
                            font.label = c("bold", 10),
                            label.rectangle = TRUE,
                            xlab = bquote(~"Log"[2]~"mean expression"),
                            ylab = bquote(~"Log"[2]~"mean expression")

  )  %>%

    ggpar(
      title = title,
      submain = contrasts,
      caption = "MA Plot - LFC apelgm",
      ggtheme = theme_pubr(),
      legend.title = "Gene Expression",
      font.title = c("bold", "brown", 16),
      font.subtitle = c("bold", "dark grey", 12),
      font.caption = c("bold.italic", "royal blue", 10),
      font.legend = c("bold", "black", 8),
      font.x = c("bold", "black", "11"),
      font.y = c("bold", "black", "11"),
      font.tickslab = c("bold", "black", "8"),
      font.family = "Calibri",

      #ylim = c(-2, 2)

    )

  #Visualise maplot


  return(maplot_apeglm)
  #ggsave(filename = paste0(plots_dir, title, "MA_Plots.png"))
}

#PCA functions ####

PCA_BIPLOT <- function(vsd, gene_prop = 1, design, batch = NA, title,
                       remove_batch_effect = FALSE){

  ##Using PCAtools

  ###Number of genes to use in PCA
  if(gene_prop < 0 | gene_prop > 1){stop("gene_prop must be between 0 and 1")}

  selectgenes <- floor(gene_prop * nrow(vsd))



  ##To correct for batch effect:
  if(is.na(batch) & remove_batch_effect){stop("Batch column must be provided to remove batch effect")}
  if(remove_batch_effect){
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd[[batch]])
  }


  ##Run PCA
  pca_data <- PCAtools::pca(mat = assay(vsd),
                            metadata = colData(vsd),
                            removeVar = (1-gene_prop))

  ##plot PC1 vs PC2
  Biplot <- biplot(pca_data,
                   showLoadings = FALSE,
                   lab = NULL,
                   #lab = pca_data$metadata$rowid,
                   #drawConnectors = TRUE,
                   #boxedLabels = TRUE,
                   #labhjust = 2,
                   #labvjust = 4,
                   colby = design,
                   #shape = NULL,
                   shape = if(!is.na(batch)){batch},
                   #encircle = TRUE,
                   #ellipse = TRUE,
                   #ellipseConf = 0.95,
                   #ellipseFill = TRUE,
                   #ellipseAlpha = 1/4,
                   #ellipseLineSize = 1.0,
                   #showLoadingsNames = TRUE,
                   #boxedLoadingsNames = TRUE,
                   #ntopLoadings = 10,
                   #sizeLoadingsNames = 2,
                   #colLoadingsNames = "black",
                   #fillBoxedLoadings = "white",
                   #drawConnectorsLoadings = TRUE,
                   #alphaLoadingsArrow = 1,
                   #lengthLoadingsArrowsFactor = 1,
                   #legendPosition = 'right',
                   legendLabSize = 10,
                   legendIconSize = 4.0) %>%

    ggpar( title = title,
           subtitle = "PCA Plot - all samples",
           caption = paste0("PCA using: ", selectgenes, " genes") ,
           legend = "right",
           legend.title = list(color = design, shape = batch),
           palette = c("red4", "blue4", "green4", "orange", "black", "brown", "violet"),
           shapekey = c(20, 4, 6),
           font.title = c("bold", "brown", 18),
           font.subtitle = c("bold", "dark grey", 12),
           font.caption = c("bold.italic", "royal blue", 10),
           font.legend = c("bold", "black", 8),
           font.x = c("bold", "black", "11"),
           font.y = c("bold", "black", "11"),
           font.tickslab = c("bold", "black", "8"),
           font.family = "Calibri",
           x.text.angle = 0,
           ggtheme = theme_pubr()
    )
  #Bp <<- Biplot
  #return(Bp)
  ifelse(remove_batch_effect,
         Biplot$labels$subtitle <- paste0(Biplot$labels$subtitle, ", batch corrected"),
           FALSE)
  return(Biplot)
  #ggsave(paste0(plots_dir, title, ifelse(remove_batch_effect, "_batch_corrected_", "_"), "PCA_all_samples.png"))
}

Plot_pca_loadings <- function(vsd,
                              gene_prop = 1,
                              nPCs = 2,
                              title = Plot_title,
                              batch = NA,
                              remove_batch_effect = FALSE
                              ) {

  ##Using PCAtools

  ###Number of genes to use in PCA
  if(gene_prop < 0 | gene_prop > 1){stop("gene_prop must be between 0 and 1")}

  selectgenes <- floor(gene_prop * nrow(vsd))

  ##To correct for batch effect:
  if(is.na(batch) & remove_batch_effect){stop("Batch column must be provided to remove batch effect")}
  if(remove_batch_effect){
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd[[batch]])
  }




  ##Run PCA
  pca_data <- PCAtools::pca(mat = assay(vsd),
                            metadata = colData(vsd),
                            removeVar = (1-gene_prop))

  ##Create a data.frame of selected PC loadings
  Loads_PCs <- PCAtools::getLoadings(pca_data, 1) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl_gene_id")

  for (i in 2:nPCs) {
    Loads_PCs <- Loads_PCs %>% dplyr::mutate(getLoadings(pca_data, i))

  }

  ##save in results folder
  write.table(Loads_PCs,
              paste0(results_dir, title, ifelse(remove_batch_effect, "_batch_corrected_", "_"),
                     "PCA_loadings.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              row.names=FALSE,
              quote=FALSE)

  ##Make a loadings plot for selected loadings
  PLoadings <- plotloadings(pca_data,
                            components = getComponents(pca_data, 1:nPCs),
                            rangeRetain = 0.1,
                            labSize = 3.0,
                            #title = 'Loadings plot',
                            #subtitle = 'PC1, PC2',
                            #caption = 'Top 10% variables',
                            shape = 23,
                            shapeSizeRange = 5,5,
                            col = c("dark red", "violet", "dark blue"),
                            colMidpoint = 0,

                            drawConnectors = TRUE,
                            lengthConnectors = unit(0.005, "npc"),
                            positionConnectors = "right",
                            #labvjust = 0.5,

  ) %>%


    ggpar( title = title,
           subtitle = "All samples",
           caption = "PC Loadings - Top 10% variables: ",
           legend = "right",
           font.title = c("bold", "brown", 18),
           font.subtitle = c("bold", "dark grey", 14),
           font.caption = c("bold.italic", "royal blue", 10),
           font.legend = c("bold", "black", 8),
           font.x = c("bold", "black", "11"),
           font.y = c("bold", "black", "11"),
           font.tickslab = c("bold", "black", "8"),
           font.family = "Comic Sans MS",
           x.text.angle = 0,



           ggtheme = theme_pubr()

    )


  ##Re-label the loadings with gene symbols:
  ##read in the annotated results table
  results_annotated <- read.delim(paste0(results_dir,"deseq2_results_annotated.tsv")) %>%
    as.data.frame()

  #pull the loadings labels and join with annotated table
  a <- dplyr::select(results_annotated, ensembl_gene_id, hgnc_symbol)
  b <- dplyr::select(PLoadings$data, var) %>%
    left_join(a,
              by=c("var" = "ensembl_gene_id"))

  ##edit loadings label to be symbols
  PLoadings$layers[[3]]$mapping$label <- as.character(b$hgnc_symbol)
  PLoadings$layers[[3]]$geom$default_aes$fontface = "bold"

  ifelse(remove_batch_effect,
         Biplot$labels$subtitle <- paste0(Biplot$labels$subtitle, ", batch corrected"),
         FALSE)


  return(PLoadings)
  ##ggsave(paste0(plots_dir, title, ifelse(remove_batch_effect, "_batch_corrected", FALSE), "_PCA_all_samples.png"))

}


##Volcano Plots ####

Volcano_Plots <- function(de_seq,
                          test,
                          res,
                          title = Plot_title,
                          labels = FALSE,
                          pval = 1e-6,
                          FC = 1.5,
                          color_specifc = NULL,
                          label_list = NULL,
                          filter_specific = NULL){

  contrast_list <- DESeq2::resultsNames(de_seq)
  contrasts <- contrast_list[str_which(contrast_list, test)]
  #res <- DESeq2::results(de_seq, name = contrasts, tidy = FALSE, independentFiltering = TRUE)

  ##read in annotated results table
  results_annotated <- read.delim(paste0(results_dir,"deseq2_results_annotated.tsv")) %>%
    as.data.frame()


  resLFC <- lfcShrink(de_seq, coef= contrasts, type="apeglm", res = res) %>%
    as.data.frame() %>%
    #rownames_to_column("ensembl_gene_id") %>%
    dplyr::mutate(hgnc_symbol = results_annotated$hgnc_symbol,
                  ENTREZID = results_annotated$ENTREZID)

  #pathways
  ##HALLMARKS GENE SETS
  HS_HALLMARK <- msigdbr(species = "Homo sapiens", category = "H")
  HS_HALLMARK <- split(x = HS_HALLMARK$gene_symbol, f = HS_HALLMARK$gs_name)


  #filter specific

  #resLFC <- resLFC %>% filter(hgnc_symbol %in% HS_HALLMARK$HALLMARK_TNFA_SIGNALING_VIA_NFKB)

  if (!is.null(filter_specific)) {
    resLFC <- resLFC %>% filter(hgnc_symbol %in% filter_specific)
  }



  #Defining colour values:
  if (!is.null(color_specifc)) {
    color_specific <- data.frame(color = color_specific) %>% rownames_to_column("hgnc_symbol")
    resLFC <- resLFC %>% left_join(color_specific) %>% mutate(color = tidyr::replace_na(color, "black"))
    keyvals <- resLFC$color
    names(keyvals) <- resLFC$color
    names(keyvals)[keyvals == 'black'] <- 'no module'
  } else {
  Low_colour <- resLFC$log2FoldChange < -FC & resLFC$padj < pval

  High_colour <- resLFC$log2FoldChange > FC & resLFC$padj < pval

  keyvals <- ifelse(
    High_colour, 'dark green',
    ifelse(Low_colour, 'red',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'dark green'] <- 'high'
  names(keyvals)[keyvals == 'black'] <- 'mid'
  names(keyvals)[keyvals == 'red'] <- 'low' }


  #Defining which genes to label


  if (!is.null(label_list)) {
    label_LFC <- resLFC %>% filter(hgnc_symbol %in% label_list)
  } else
    {label_LFC <- results_annotated %>%
    dplyr::filter(abs(log2FoldChange) > FC, padj < pval)}

  ifelse(labels,
         labs <-label_LFC$hgnc_symbol,
         labs <- ""

         )


  Volcano <- EnhancedVolcano(resLFC,
                             lab = resLFC[["hgnc_symbol"]],
                             #lab = rownames(resLFC),
                             selectLab = labs,
                             #selectLab = ifelse(labels, resLFC$hgnc_symbol[which(names(keyvals) %in% c("high", "low"))], ""),
                             title = title,
                             subtitle = paste0(contrasts, "_using_FC_", FC, "_padj_", pval),
                             labSize = 2.0,
                             labFace = "bold.italic",
                             labCol = "dark blue",
                             x = 'log2FoldChange',
                             y = 'padj',
                             xlab = bquote(~Log[2]~ 'fold change'),
                             ylab = bquote(~-Log[10]~ 'p-value change'),
                             xlim = c(-10, 10),
                             ylim = c(0, 45),
                             pCutoff = pval,
                             FCcutoff = FC,
                             #pointSize = 2.0,
                             colCustom = keyvals,
                             pointSize = ifelse(rownames(resLFC) %in% label_LFC[, 1], 2, 1),
                             #pointSize = c(ifelse(abs(res$log2FoldChange)>5, 3, 1)),
                             gridlines.major = FALSE,
                             gridlines.minor = FALSE,
                             colAlpha = 4/5,
                             legendPosition = 'right',
                             legendLabSize = 6,
                             legendIconSize = 2.0,
                             drawConnectors = TRUE,
                             widthConnectors = 0.8,
                             boxedLabels = TRUE


  )

  #options(ggrepel.max.overlaps = Inf)

  Volcano <- ggpar(Volcano,
                   legend.title = "Differentially expressed genes",
                   ggtheme = theme_pubclean(),
                   caption = paste("Volcano plot using", Volcano$labels$caption, sep = " "),
                   font.caption = c("bold.italic", "royal blue", 12),
                   font.legend = c(8, "bold"),
                   font.main = c(24, "bold", "dark red"),
                   font.submain = c(14, "bold.italic"))





  return(Volcano)



}


### GSEA ####


GSEA_plots <-  function(pathways,
                        res = paste0(results_dir, "deseq2_results_annotated.tsv"),
                        title,
                        control,
                        test){

  caption <- str_to_upper(pathways)
  ##read in the annotated results table - should have the same column headings as DeSeq2 results table
  if(is.character(res)){
    results_annotated <- read.delim(res) %>%
      as.data.frame()
  }

  if(class(res) == "DESeqResults"){
    results_annotated <- as.data.frame(res) %>%
      tibble::rownames_to_column("ensembl_gene_id")

  }




  ##read in the annotated results table
 # results_annotated <- read.delim(paste0(results_dir,"deseq2_results_annotated.tsv")) %>%
  #  as.data.frame()

  ##create ranked gene list for fgsea
  gseaInput <-  results_annotated %>%
    dplyr::filter(!is.na(hgnc_symbol), !is.na(stat)) %>%
    arrange(stat)

  ranks <- pull(gseaInput,stat)
  names(ranks) <- gseaInput$hgnc_symbol


  ###choose pathway databases to load.

  ##REACTOME GENE SETS
  HS_CP_REACTOME <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  HS_CP_REACTOME <- split(x = HS_CP_REACTOME$gene_symbol, f = HS_CP_REACTOME$gs_name)


  ##HALLMARKS GENE SETS
  HS_HALLMARK <- msigdbr(species = "Homo sapiens", category = "H")
  HS_HALLMARK <- split(x = HS_HALLMARK$gene_symbol, f = HS_HALLMARK$gs_name)

  ##KEGG GENE SETS
  HS_CP_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  HS_CP_KEGG <- split(x = HS_CP_KEGG$gene_symbol, f = HS_CP_KEGG$gs_name)


  ##Choose gene sets list
  #pathways <- HS_HALLMARK
  c <- str_to_upper(pathways)
  pathway_dbs <- list("HS_HALLMARK" = HS_HALLMARK, "HS_CP_KEGG" = HS_CP_KEGG, "HS_CP_REACTOME" = HS_CP_REACTOME)
  a <-c("HS_HALLMARK", "HS_CP_KEGG", "HS_CP_REACTOME")
  b <- a[str_which(a, str_to_upper(pathways))]
  pathways <- pathway_dbs[[b]]
  ##GSEA analysis using the "ranks" table against the pathways in the database.

  fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)



  ## Show results in a nice table

  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::mutate(Pathway = str_replace_all(pathway, "_", " "),
                  Adjusted_p_value = padj,
                  Normalised_Enrichment_Score = NES,
                  Size = size,
                  Key_genes = leadingEdge) %>%
    dplyr::mutate(Pathway = str_replace(Pathway, c, "")) %>%
    dplyr::select(Pathway, Normalised_Enrichment_Score, Size, Adjusted_p_value)
  #%>% dplyr::filter(abs(Normalised_Enrichment_Score) > 2.3)

  ##save results table
  #return(fgseaResTidy)
  write.table(fgseaResTidy,
              file = paste0(results_dir, "fgsea_", b ,"_pathways.tsv"),
              col.name=TRUE,
              sep="\t",
              na = "NA",
              append = FALSE,
              row.names=FALSE,
              quote=FALSE)




  ##kable table of results
  FGSEA_Results <- fgseaResTidy %>% kbl(col.names = c("Pathway", "Normalised Enrichment Score", "Size", "Adjusted p-value"),
                                        caption = paste0(Plot_title, " ", Control, "_vs_", Test, " Top Differentially Expressed Pathways"),
                                        digits = 2,
                                        escape = F,
                                        align = "c"
  ) %>%
    kable_paper(html_font = "Comic Sans MS",
                full_width = T, bootstrap_options = c("striped", "bordered", "hover")) %>%
    column_spec(2, color = "white", background = ifelse(fgseaResTidy$Normalised_Enrichment_Score < 0, "red", "green")) %>% footnote("Using Kable")



  ###Plots enrichment scores

  FGSEA_plot <- (ggbarplot(fgseaResTidy %>% dplyr::arrange(abs(Normalised_Enrichment_Score)) %>% head(25) %>% dplyr::arrange(desc(Normalised_Enrichment_Score)),
                           x = "Pathway",
                           y = "Normalised_Enrichment_Score",
                           fill = "Adjusted_p_value",
                           #fill = "Size",
                           #palette = "jco",
                           sort.val = "desc",
                           label = FALSE,
                           xlab = FALSE,
                           size = 1,
                           width = 1,
  )
  + scale_x_discrete(labels = function(x) str_wrap(x, 15))

  + scale_fill_steps(
    low = "#132B43",
    high = "#56B1F7",
    space = "Lab",
    na.value = "grey50",
    guide = "coloursteps",
    aesthetics = "fill",
    name = "P-value",
    n.breaks = 5
  )) %>%

    ggpar(#legend.title = list(color = "P-values"),
      rotate = TRUE,
      title = title,
      submain = paste0(control, " vs ", test),
      caption = paste0("GSEA with ", caption, " gene sets from MSigDB"),
      ggtheme = theme_pubr(),
      #ggtheme = theme(legend.direction = "vertical"),
      xlab = "",
      ylab = "Normalized Enrichment Score",
      font.title = c("bold", "brown", 18),
      font.subtitle = c("bold", "dark grey", 14),
      font.caption = c("bold.italic", "royal blue", 12),
      font.legend = c("bold", "black", 6),
      font.x = c("bold", "black", "11"),
      #x.text.angle = 90,
      font.y = c("bold", "black", "11"),
      font.tickslab = c("bold", "black", "8"),
      font.family = "Calibri",

    )

  FGSEA_plot
  ggsave(filename =  paste0(plots_dir, Plot_title, "_", c, "_FGSEA_plot.png"), width = 10, height = 14, units = "in")
  GSEA_return <- list(FGSEA_Results, FGSEA_plot)
  return(GSEA_return)

}


##ORA analysis ####

ORA_GoSeq <- function(res = paste0(results_dir, "deseq2_results_annotated.tsv"),
                      pvalue = 0.05,
                      FC = 1,
                      title,
                      control,
                      test){

  ##read in the annotated results table - should have the same column headings as DeSeq2 results table
  if(is.character(res)){
  results_annotated <- read.delim(res) %>%
    as.data.frame()
  }

  if(class(res) == "DESeqResults"){
    results_annotated <- as.data.frame(res) %>%
      tibble::rownames_to_column("ensembl_gene_id")

  }
    #First, create list of DEGs with chosen pvalue and fold change cut-off

  if (FC > 0) {genes <- results_annotated$padj < pvalue &
    results_annotated$log2FoldChange > FC &
    !is.na(results_annotated$padj)
  }

  if (FC < 0) {genes <- results_annotated$padj < pvalue &
    results_annotated$log2FoldChange < FC &
    !is.na(results_annotated$padj)
  }


  names(genes) <- results_annotated[, 1]


  ##Then run Probability Weighting Function using hg38 dataset - TxDb.Hsapiens.UCSC.hg38.knownGene - this database is used to calculate the transcript length - goseq corrects transcript length using gene lengths from the database.
  ## need to use ensGene as gene id here as ENTREZ "geneid" and symbols "geneSymbol" have duplicates.

  pwf <- nullp(genes, "hg38", "ensGene")

  goseq_res <- goseq::goseq(pwf,
                     "hg38",
                     "ensGene",
                     #test.cats = c("GO:BP")
                     test.cats= c("GO:CC", "GO:BP", "GO:MF", "KEGG")
                     )



  ##Plot the top 25 results pvalue:
  goPlot <- ((goseq_res %>%
               top_n(25, wt=-over_represented_pvalue) %>%
               mutate(hitsPerc=numDEInCat*100/numInCat) %>%
               slice_max(hitsPerc, prop = 0.5) %>%
               ggplot(aes(x=hitsPerc,
                          y=term,
                          colour=over_represented_pvalue,
                          size=numDEInCat)) +
               geom_point() +
               expand_limits(x=0) +
               labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

 )  +   scale_y_discrete(labels = function(x) str_wrap(x, 20))) %>%

    ggpar(legend.title = list(color = "P-values", size = "Count"),
          rotate = FALSE,
          title = title,
          submain = paste0(control, " vs ", test),
          caption = "Go seq plot",
          ggtheme = theme_pubr(),
          legend = "right",
          xlab = "Hits (%)",

          ylab = "GO term",
          font.title = c("bold", "brown", 18),
          font.subtitle = c("bold", "dark grey", 14),
          font.caption = c("bold.italic", "royal blue", 12),
          font.legend = c("bold", "black", 8),
          font.x = c("bold", "black", "11"),
          font.y = c("bold", "black", "11"),
          font.tickslab = c("bold", "black", "8"),
          font.family = "Calibri",

    )



  ###save a list of significant DEGs for use in online tools
  universe <- results_annotated %>% dplyr::select(ensembl_gene_id)

  if (FC < 0) {sigGenes <- results_annotated %>%
    dplyr::filter(padj < pvalue, log2FoldChange < FC
                  #!is.na(ENSEMBL)
    ) %>% dplyr::select(ensembl_gene_id)
  }

  if (FC > 0) {sigGenes <- results_annotated %>%
    dplyr::filter(padj < pvalue, log2FoldChange > FC
                  #!is.na(ENSEMBL)
    ) %>% dplyr::select(ensembl_gene_id)
  }


  ##list of all genes in assay
  write.table(universe,
              paste0(results_dir, title, "_gene_universe_ensemblid.txt"),
              col.name=FALSE,
              sep="\t",
              na = "NA",
              row.names=FALSE,
              quote=FALSE)

  ##choose txt file title
  over_or_under <- ifelse(FC > 0, "over", "under")

  ##list of significant genes with chosen cut offs.
  write.table(sigGenes,
              paste0(results_dir, "gene_DEG_log2fold", over_or_under, FC,"_pvalue_", pvalue, "_ensemblid.txt"),
              col.name=FALSE,
              sep="\t",
              na = "NA",
              row.names=FALSE,
              quote=FALSE)


  ggsave(filename =  paste0(plots_dir, Plot_title, "_", over_or_under, "represented_GO", "_plot.png"), width = 10, height = 14, units = "in")



  ORA_return <- list(goseq_res, goPlot, universe, sigGenes)
  return(goPlot)
  #return(ORA_return)
}






### clusterProfiler ####

clusterProfiler_Plots <- function(res = paste0(results_dir, "deseq2_results_annotated.tsv"),
                                  title) {


  ##Analysis with clusterProfiler

  results_annotated <- read.delim(res)
  universe <- results_annotated %>% pull(ensembl_gene_id)
  sigGenes <- results_annotated %>%
    dplyr::filter(padj < 0.05, log2FoldChange < -1.5
                  #!is.na(ENSEMBL)
    ) %>% pull(ensembl_gene_id)

  enrich_go <- enrichGO(
    gene= sigGenes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    universe = universe,
    qvalueCutoff = 0.05,
    readable=TRUE
  ) %>% pairwise_termsim()


  ##ggplot(enrich_go, aes(Count)) + geom_density(bw = "SJ")
  ##Visualise using dotplot and emapplot

  dotplot <- (dotplot(enrich_go) +
                scale_colour_gradientn(limits=c(0, 0.05), colours=rev(brewer.pal(9,"YlOrRd")), name ="P-value") + scale_y_discrete(labels = function(x) str_wrap(x, width = 20))
              + scale_size(range = c(1, 6), name = "Number of genes", breaks = seq(0, max(enrich_go$Count), by=20))
  ) %>%
    ggpar(
      xlab = "Gene Ratio",
      ylab = "Pathways",
      title = Plot_title,
      submain = "Gene Set Enrichment",
      caption = "Enriched Pathways",
      #palette = inferno(10),
      ggtheme = theme_pubr(),
      legend = "top",
      font.title = c("bold", "brown", 18),
      font.subtitle = c("bold", "dark grey", 14),
      font.caption = c("bold.italic", "royal blue", 12),
      font.legend = c("bold", "black", 6),
      font.x = c("bold", "black", "11"),
      font.y = c("bold", "black", "11"),
      font.tickslab = c("bold", "black", "6"),
      font.family = "Comic Sans MS",
      x.text.angle = 0
    )




  #enrich_go2 <- pairwise_termsim(enrich_go)
  emapplot <- (emapplot(enrich_go,
                        showCategory = 20,
                        color = "p.adjust",
                        layout = "nicely",
                        #legend_n = 5,
                        cex_label_category = 0.5,
                        cex_category = 2,

  )

  + scale_colour_gradientn(limits=c(0, 0.05), colours=rev(brewer.pal(9,"Blues")), name ="P-value", )
  + scale_size(range = c(2, 8), name = "Number of genes", breaks = seq(0, max(enrich_go$Count), by=20))
  ##guide = guide_colorbar(reverse = TRUE)
  + guides(colour = guide_legend(title.position = "top"))

  ##  + geom_node_label(aes(label = str_wrap(name, 15)),repel=TRUE, fontface = "bold", family = "Comic Sans MS", size = 2, color = "red", check_overlap = TRUE, segment.colour = "#666666", segment.size = 0.5, arrow = NULL, force = 10, alpha = 0.1,) ###this can be used to add labels manually, use gginnards or ggedit to remove existing layer

  ) %>%
    ggpar(

      title = title,
      submain = "Gene Set Enrichment",
      caption = "Enriched Pathways - Emapplot",
      #palette = inferno(10),
      #ggtheme = theme_grey(),
      legend = "right",
      font.title = c("bold", "brown", 18),
      font.subtitle = c("bold", "dark grey", 14),
      font.caption = c("bold.italic", "royal blue", 12),
      font.legend = c("bold", "black", 6),
      #font.x = c("bold", "black", "11"),
      #font.y = c("bold", "black", "11"),
      #font.tickslab = c("bold", "black", "6"),
      tickslab = FALSE,
      #font.family = "Comic Sans MS",
      x.text.angle = 0
    )

  ## resetting some pathway name displays to be readable
  emapplot$data$name <- str_wrap(emapplot$data$name, 10)
  emapplot$layers[[3]]$geom$default_aes$fontface <- "bold"



  gene_set_plots <- ggarrange(dotplot, emapplot, nrow = 2, ncol = 1)
  return(gene_set_plots)
  #dotplot
  ##Alternatively use the genes list to write a text file using write.table to create to files - a list of all genes in the analysis and a list of significant differentially expressred genes.These text files can be analysed in GORilla



  }





### make heatmap ####
make_heatmap <- function(res = paste0(results_dir, "deseq2_results_annotated.tsv"),
                         pathways,
                         vsd = DeSeq2_return[["vsd"]],
                         title = Plot_title
) {

  results_annotated <- read.delim(res)


  ##Analyze one of the pathways.


  ###choose pathway databases to load.

  ##REACTOME GENE SETS
  HS_CP_REACTOME <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  HS_CP_REACTOME <- split(x = HS_CP_REACTOME$gene_symbol, f = HS_CP_REACTOME$gs_name)


  ##HALLMARKS GENE SETS
  HS_HALLMARK <- msigdbr(species = "Homo sapiens", category = "H")
  HS_HALLMARK <- split(x = HS_HALLMARK$gene_symbol, f = HS_HALLMARK$gs_name)

  ##KEGG GENE SETS
  HS_CP_KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  HS_CP_KEGG <- split(x = HS_CP_KEGG$gene_symbol, f = HS_CP_KEGG$gs_name)


  ##Choose gene sets list
  #pathways <- HS_HALLMARK
  pathway_dbs <- list("HS_HALLMARK" = HS_HALLMARK, "HS_CP_KEGG" = HS_CP_KEGG, "HS_CP_REACTOME" = HS_CP_REACTOME)
  a <-c("HS_HALLMARK", "HS_CP_KEGG", "HS_CP_REACTOME")
  b <- a[str_which(a, str_to_upper(pathways))]
  pathways <- pathway_dbs[[b]]



  ##First create a list of genes by identifying all the genes in the de-seq results table that are in the chosen Pathway.

  my_genes <- dplyr::filter(results_annotated, hgnc_symbol %in% pathways[["REACTOME_TRANSLATION"]]) %>%
    pull(ensembl_gene_id)

  #my_genes <- head(my_genes, 50)


  ###pull of the log normalised counts of selected genes into a matrix.
  mat <- assay(vsd)[my_genes,]

  rownames(mat) <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = rownames(mat),
    keytype = "GENEID",
    column = "SYMBOL"
  )

  ###Subtract each sample's gene count from the mean count for that genes across all sames. Now "mat" contains deviation values from the mean expression of each gene.

  mat <- mat - rowMeans(mat)


  ##read in the metadata data.frame.
  sampleinfo_for_gsea <-
    as.data.frame(colData(vsd))[c(deseq2_design_condition)]

  ##Rename column names from samplenames to IDs.
  colnames(mat) <- rownames(sampleinfo_for_gsea) <- vsd$ID

  ###Plot Heatmap to visualise clustering of samples as per differential expression of this gene set.



  #quantile_breaks <- function(xs, n = 10) {
  #  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  #  breaks[!duplicated(breaks)]
  #}

  #mat_breaks <- quantile_breaks(mat, n = 11)




  my_colour <- list(Samples = c("#5977ff", "#f74747"))
  names(my_colour$Samples) <- c(Control, Test)

    #Batch = c("1" = "#82ed82", "2" = "#9e82ed"),
    #gleason_Score = c("G6" = "#e89829", "G9" = "#cc4ee0", "Normal" = "royal blue")




  PMAP <- pheatmap(mat,
                   annotation_col = sampleinfo_for_gsea,
                   annotation_colors = my_colour,
                   #main = "Default Heatmap",
                   color = colorRampPalette(c("dark red", "white", "dark green"))(2),
                   #color = inferno(length(mat_breaks) - 1),
                   #breaks = seq(min(mat), max(mat), length.out = 3),
                   breaks = c(min(mat), 0, max(mat)),
                   cellheight = 5,
                   show_colnames     = TRUE,
                   show_rownames     = TRUE,
                   fontsize_row = 5,
                   fontsize_col = 8,

  ) %>% ggplotify::as.ggplot() %>%

    ggpar( ticks = FALSE,
           tickslab = FALSE,
           xlab = FALSE,
           ylab = FALSE,
           title = title,
           subtitle = paste(Control, "vs", Test),
           caption = b,
           font.title = c("bold", "brown", 18),
           font.subtitle = c("bold", "dark grey", 14),
           font.caption = c("bold.italic", "royal blue", 10),

           #font.family = "Comic Sans MS",


    )






  return(PMAP)


}

#### makde module heatmap for wgcna ####

make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.

  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("refinebio_accession_code")

  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(refinebio_accession_code, time_point, refinebio_subject) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "refinebio_accession_code") %>%
    # Arrange by patient and time point
    dplyr::arrange(time_point, refinebio_subject) %>%
    # Store sample
    tibble::column_to_rownames("refinebio_accession_code")

  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    time_point = col_annot_df$time_point,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(time_point = c("recovering" = "#f1a340", "acute illness" = "#998ec3"))
  )

  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)

  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()

  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()

  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )

  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )

  # Return heatmap
  return(heatmap)

}


#### ####


gnet2 <- function(input,reg_names,init_method= 'boosting',init_group_num = 4,max_depth = 3,
                 cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5,
                 heuristic = TRUE,max_group = 0,force_split = 0.5,nthread = 4){
  if(is(input,class2 = "SummarizedExperiment")){
    input <- assay(input)
  }
  if(max_group == 0){
    max_group <- init_group_num
  }
  input <- input[apply(input, 1, var) > 0.0001,]
  gene_data <- input[!rownames(input)%in%reg_names,,drop=FALSE]
  regulator_data <- input[rownames(input)%in%reg_names,,drop=FALSE]
  result_all <- GNET2:::run_gnet(gene_data,regulator_data,init_method,init_group_num,max_depth,cor_cutoff,
                         min_divide_size,min_group_size,max_iter,heuristic,max_group,force_split,nthread = nthread)
  reg_group_table <- result_all[[1]]
  gene_group_table <- result_all[[2]]

  # sanity check: remove all modules without any genes assigned



  group_idx <- unique(reg_group_table[,1])
  reg_group_table_filtered <- gene_group_table_filtered <- NULL
  avg_cor_list <- c()
  current_group_idx <- 1
  regulators <- target_genes <- list()
  for (i in seq_len(length(group_idx))) {
    if(sum(gene_group_table[,2]==group_idx[i])>=min_group_size){
      current_tree <- reg_group_table[reg_group_table[,1] == group_idx[i],]
      current_tree[,1] <- current_group_idx
      current_gene_group <- gene_group_table[gene_group_table$group==group_idx[i],]
      current_gene_group$group<- current_group_idx
      reg_group_table_filtered <- rbind(reg_group_table_filtered,current_tree)
      gene_group_table_filtered <- rbind(gene_group_table_filtered,current_gene_group)

      cor_m <- cor(t(gene_data[current_gene_group$gene,,drop=FALSE]))
      avg_cor_list <- c(avg_cor_list,mean(cor_m[upper.tri(cor_m)]))

      regulators[[i]] <- rownames(regulator_data)[current_tree[,2]+1]
      target_genes[[i]] <- current_gene_group$gene
      current_group_idx <- current_group_idx + 1
    }
  }

  if(nrow(reg_group_table_filtered)<=2)warning('Too few modules generated, you may wish to try with higher cor_cutoff.')
  return(list('gene_data' = gene_data,'regulator_data' = regulator_data,'group_score' = avg_cor_list,
              'reg_group_table' = reg_group_table_filtered,'gene_group_table' = gene_group_table_filtered,
              'modules_count' = current_group_idx-1,'regulators' = regulators,
              'target_genes' = target_genes))
}
