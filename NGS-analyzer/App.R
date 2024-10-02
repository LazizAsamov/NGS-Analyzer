# Clear the environment
rm(list=ls())

# Load required libraries
library(shiny)
library(dplyr)
library(DESeq2)
library(shinyjs)
library(future)
library(promises)
library(shinyjqui)
library(bslib)
library(DT)
library(EnhancedVolcano)
library(ggplot2)
library(plotly)
library(Glimma)
library(pheatmap)

# Set the working directory
#setwd("")

# Define UI
ui <- htmlTemplate(
  filename = "www/custom_layout.html",
  
  # Include ShinyJS for JavaScript operations
  ShinyJS = useShinyjs(),
  
  # Inputs for specifying conditions and uploading files
  NumInput = numericInput("numConditions", "How many conditions do you have?", value = 2, min = 1),
  ConditionUI = uiOutput("conditionInputs"),
  FileUI = uiOutput("fileInputs"),
  
  # Results and plots
  BTAB = DTOutput('dds_result', width = "80%"),
  p_his= plotlyOutput('prob_his', width = "100%", height = "100%"),
  disp_plt= plotlyOutput('disp_plt', width = "100%", height = "100%"),
  plt_volcano = plotOutput('volcano', width = "100%", height = "100%"),
  pca_plt = plotlyOutput('pca_plt', width = "100%", height = "100%"),
  GeneSelect = uiOutput('GeneSelect'),
  count_plt = plotlyOutput('count_plt', width = "100%", height = "80%"),
  heat_plt = plotlyOutput('heatmap', width = "100%", height = "100%"),
  
  # Sliders for threshold settings
  SliderPvalue = sliderInput("sliderP", label = "Set p-value threshold", min = 0, max = 1, value = 0.05, step= 0.01),
  SliderLog2F = sliderInput('sliderLog', label = "Set Log2FC threshold", min = 0, max = 3, value = 3)
)

# Define Server
server <- function(input, output, session) {
  
  # Enable sortable elements
  jqui_sortable("#sortable_container")
  
  # Set maximum upload size
  options(shiny.maxRequestSize = 8000 * 1024^2)
  
  # Reactive values to store conditions, files, and metadata
  conditions_list <- reactiveVal()
  user_files <- reactiveVal(list()) 
  metadata_df <- reactiveVal(data.frame(condition = character(), row.names = NULL, stringsAsFactors = FALSE))
  
  # Dynamically generate condition inputs based on user input
  output$conditionInputs <- renderUI({
    req(input$numConditions)
    tagList(
      lapply(1:input$numConditions, function(i) {
        textInput(paste0("condition_", i), label = paste("Condition", i), placeholder = "e.g., Treated, Untreated")
      }),
      actionButton("submitConditions", "Submit Conditions", class = "btn btn-light")
    )
  })
  
  # Generate file inputs after user specifies conditions
  observeEvent(input$submitConditions, {
    # Collect conditions entered by the user
    conditions <- sapply(1:input$numConditions, function(i) input[[paste0("condition_", i)]])
    conditions_list(conditions)
    
    # Render file input for each condition
    output$fileInputs <- renderUI({
      req(conditions_list())
      
      tagList(
        lapply(conditions, function(condition) {
          fileInput(paste0("file_", condition), label = paste("Upload files for", condition), multiple = TRUE)
        }),
        fileInput("genome", label = 'Upload reference genome'),
        fileInput("annotation", label = 'Upload annotation file'),
        actionButton("submitFiles", "Submit Files", class= 'btn btn-light'),
        actionButton("gotoResults", "Go to Results", class = "btn btn-light"),
        actionButton("run_pipeline", "Run Pipeline", class = "btn btn-danger")
      )
    })
  })
  
  # Store metadata and upload files upon file submission
  observeEvent(input$submitFiles, {
    req(conditions_list)
    conditions <- conditions_list()
    
    # Parse file prefixes and update metadata
    for(condition in conditions) {
      files <- input[[paste0("file_", condition)]]$name
      if (is.null(files)) return(NULL)
      
      for (file_name in files) {
        if (grepl("_r1.fastq$", file_name)) {
          existing_metadata <- metadata_df()
          prefix <- sub("_r1.fastq$", "", file_name)
          new_metadata <- data.frame(condition = condition, row.names = prefix, stringsAsFactors = FALSE)
          combined_df <- rbind(existing_metadata, new_metadata)
          metadata_df(combined_df)
        }
      }
    }
    
    # Define destination directories for storing uploaded files
    dest_dir <- "./data"
    if (!dir.exists(dest_dir)) dir.create(dest_dir)
    
    # Move condition-specific files to destination folder
    for (condition in conditions) {
      file_info <- input[[paste0("file_", condition)]]
      if (!is.null(file_info)) {
        for (i in seq_along(file_info$name)) {
          tempfile <- file_info$datapath[i]
          filename <- file_info$name[i]
          dest_file <- file.path(dest_dir, filename)
          file.rename(tempfile, dest_file)
        }
      }
    }
    
    # Move genome and annotation files to their respective folder
    dest_dir_gen_annot <- "./genome"
    if (!dir.exists(dest_dir_gen_annot)) dir.create(dest_dir_gen_annot)
    
    file_info_gen <- input[["genome"]]
    if (!is.null(file_info_gen)) {
      for (i in seq_along(file_info_gen$name)) {
        tempfile <- file_info_gen$datapath[i]
        filename <- file_info_gen$name[i]
        dest_file <- file.path(dest_dir_gen_annot, filename)
        file.rename(tempfile, dest_file)
      }
    }
    
    file_info_annot <- input[["annotation"]]
    if (!is.null(file_info_annot)) {
      for (i in seq_along(file_info_annot$name)) {
        tempfile <- file_info_annot$datapath[i]
        filename <- file_info_annot$name[i]
        dest_file <- file.path(dest_dir_gen_annot, filename)
        file.rename(tempfile, dest_file)
      }
    }
  })
  
  # Render and print conditions and file metadata
  output$outputConditionsFiles <- renderPrint({
    req(input$submitFiles)
    Con <- conditions_list()
    m_df <- metadata_df()
    print(Con)
    print(m_df)
  })
  
  # Pipeline execution logic
  pipeline_results <- eventReactive(input$run_pipeline, {
    withProgress(message = 'Please stand by...', min = 0, max = 1, {
      shinyjs::disable("run_pipeline")
      
      # Identify paired-end files
      files_names <- list.files(path = './data')
      r1_files <- grep("_r1", files_names, value = TRUE)
      r2_files <- grep("_r2", files_names, value = TRUE)
      
      # Create dataframe for paired-end reads
      paired_files <- data.frame(forward = character(), reverse = character(), stringsAsFactors = FALSE)
      for (r1_file in r1_files) {
        sample_id <- gsub("_r1.*", "", r1_file)
        r2_file <- r2_files[grep(sample_id, r2_files)]
        if (length(r2_file) == 1) {
          paired_files <- rbind(paired_files, data.frame(forward = r1_file, reverse = r2_file))
        }
      }
      
      # FastQC Analysis
      setProgress(value = 0.2, message = 'Running FastQC')
      for (i in 1:nrow(paired_files)){
        sample_1_r1 <- paired_files[i,1]
        sample_1_r2 <- paired_files[i,2]
        
        fastqc_script <- c(
          "#!/bin/bash",
          "echo 'Starting FastQC!'",
          "pwd",
          sprintf("fastqc -o fastQC data/%s data/%s", sample_1_r1, sample_1_r2)
        )
        writeLines(fastqc_script, con = "./bash_scripts/run_fastqc.sh")
        system("chmod +x ./bash_scripts/run_fastqc.sh")
        system("./bash_scripts/run_fastqc.sh")
        system("rm ./fastQC/*.zip ; rm ./bash_scripts/run_fastqc.sh")
      }
      shinyjs::enable("run_pipeline")
      output$execStatus <- renderText("Calculation finished")
      
      # Indexing genome using BWA
      setProgress(value = 0.4, message = 'Indexing genome')
      bwa_index_script <- c(
        "#!/bin/bash",
        "source /Users/lazizasamov/micromamba/etc/profile.d/micromamba.sh",
        "micromamba activate RNA-seq",
        "export PATH=$PATH:/Users/lazizasamov/micromamba/envs/RNA-seq/bin",
        "echo 'Starting indexing!'",
        "bwa index -p indexed_genome/indexed genome/*.fna"
      )
      writeLines(bwa_index_script, con = "./bash_scripts/run_bwa_index.sh")
      system("chmod +x ./bash_scripts/run_bwa_index.sh")
      system("./bash_scripts/run_bwa_index.sh")
      system("rm ./bash_scripts/run_bwa_index.sh")
      
      # Mapping reads to reference genome
      setProgress(value = 0.6, message = 'Mapping reads to reference genome')
      for (i in 1:nrow(paired_files)){
        sample_1_r1 <- paired_files[i,1]
        sample_1_r2 <- paired_files[i,2]
        sample_prefix <- sub("_r1.fastq", "", sample_1_r1)
        
        bwa_mem_script <- c(
          "#!/bin/bash",
          "source /Users/lazizasamov/micromamba/etc/profile.d/micromamba.sh",
          "micromamba activate RNA-seq",
          "export PATH=$PATH:/Users/lazizasamov/micromamba/envs/RNA-seq/bin",
          sprintf("bwa mem indexed_genome/indexed data/%s data/%s > mapped/%s.sam", sample_1_r1, sample_1_r2, sample_prefix),
          sprintf("samtools sort -o mapped/bam/%s.bam mapped/%s.sam", sample_prefix, sample_prefix)
        )
        writeLines(bwa_mem_script, con = "./bash_scripts/run_bwa_mem.sh")
        system("chmod +x ./bash_scripts/run_bwa_mem.sh")
        system("./bash_scripts/run_bwa_mem.sh")
        system("rm ./bash_scripts/run_bwa_mem.sh")
      }
      
      # Counting reads using featureCounts
      setProgress(value = 0.8, message = 'Counting reads')
      featurecounts_script <- c(
        "#!/bin/bash",
        "source /Users/lazizasamov/micromamba/etc/profile.d/micromamba.sh",
        "micromamba activate RNA-seq",
        "export PATH=$PATH:/Users/lazizasamov/micromamba/envs/RNA-seq/bin",
        "featureCounts -s2 -p -a genome/annotationRef.gtf -o counts/count.txt -t gene mapped/bam/*.bam"
      )
      writeLines(featurecounts_script, con = "./bash_scripts/run_fc.sh")
      system("chmod +x ./bash_scripts/run_fc.sh")
      system("./bash_scripts/run_fc.sh")
      system("rm ./bash_scripts/run_fc.sh")
    })
  })
  
  # Run pipeline on button click
  observeEvent(input$run_pipeline, {
    pipeline_results()
  })
  
  # Go to results and perform differential analysis
  observeEvent(input$gotoResults, {
    # Read count data and metadata
    col_names <- c("gene_id", "chr", 'start', 'end', 'strand', 'length')
    bam_file_names <- list.files('./mapped/bam')
    for (i in bam_file_names) {
      prefix <- sub(".bam", "", i)
      col_names <- c(col_names, prefix)
    }
    count_data <- read.table("counts/count.txt", col.names = col_names)
    count_info <- read.csv("/Users/lazizasamov/Desktop/RNA-SEQ_APP/App/counts/meta_data.csv", row.names = 1)
    
    # Filter and transform data
    count_data <- count_data[-1,]
    rownames(count_data) <- NULL
    deseq_data <- dplyr::select(count_data, -c(chr, start, end, strand, length))
    rownames(deseq_data) <- deseq_data$gene_id
    deseq_data <- dplyr::select(deseq_data, -c(gene_id))
    
    # Reorder columns if needed
    if (!all(colnames(deseq_data) == rownames(count_info))) {
      deseq_data <- deseq_data[, rownames(count_info)]
    }
    
    # Convert data to numeric
    deseq_data[] <- sapply(deseq_data[], function(x) as.numeric(as.character(x)))
    
    # Create DESeq dataset
    dds <- DESeqDataSetFromMatrix(countData = deseq_data, colData = count_info, design = ~ condition)
    dds$condition <- relevel(dds$condition, ref = "untreated")
    smallestGroupSize <- 3
    keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    dds <- dds[keep,]
    
    # Run DESeq analysis
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds))
    
    # Render results table with filters
    output$dds_result <- DT::renderDT({
      tb <- subset(res, padj < input$sliderP & abs(log2FoldChange) >= input$sliderLog)
      tb <- tb %>% signif(3)
      datatable(tb, width = "10%", height = "100px")
    })
    
    # Histogram of adjusted p-values
    output$prob_his <- renderPlotly({
      his <- ggplot(res, aes(x = padj)) +
        geom_histogram(color = "white", fill = "deepskyblue3", bins = 30) +
        theme_bw() +
        labs(title = "Histogram of Adjusted p-values", x = "Adjusted p-values (padj)", y = "Counts") +
        theme(plot.title = element_text(hjust = 0.5, face="bold", size = 16))
      ggplotly(his)
    })
    
    # Dispersion plot
    output$disp_plt <- renderPlotly({
      plotDispEsts2 <- function(dds) {
        df <- as.data.frame(mcols(dds)) %>%
          dplyr::select(baseMean, dispGeneEst, dispFit, dispersion) %>%
          melt(id.vars="baseMean") %>%
          filter(baseMean > 0)
        graph <- ggplot(df, aes(x = baseMean, y = value, colour = variable)) +
          geom_point(size = 0.5) +
          scale_x_log10() +
          scale_y_log10() +
          theme_bw() +
          labs(title = "Dispersion vs BaseMean", x = "Base Mean (log scale)", y = "Dispersion (log scale)") +
          scale_colour_manual(values = c("Black", "#e41a1c", "#377eb8"), breaks = c("dispGeneEst", "dispFit", "dispersion"), labels = c("Estimate", "Fit", "Final"))
        ggplotly(graph)
      }
      plotDispEsts2(dds)
    })
    
    # Volcano plot
    output$volcano <- renderPlot({
      EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', xlim = c(min(res[["log2FoldChange"]], na.rm = TRUE), max(res[["log2FoldChange"]], na.rm = TRUE)))
    })
    
    # PCA plot
    output$pca_plt <- renderPlotly({
      vsd <- vst(dds, blind = FALSE)
      pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      pcaPlt <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed() +
        theme_minimal() +
        labs(title = "Dispersion vs BaseMean")
      ggplotly(pcaPlt)
    })
    
    # Heatmap of top genes
    output$heatmap <- renderPlotly({
      norm_counts <- counts(dds, normalized = TRUE)
      top_genes <- order(rowVars(norm_counts), decreasing = TRUE)[1:50]  # Use top 50 genes
      subset_transformed_counts <- assay(vst(dds, blind = FALSE))[top_genes, ]
      
      # Generate heatmap
      hm_plt <- plot_ly(
        x = colnames(subset_transformed_counts),
        y = rownames(subset_transformed_counts),
        z = as.matrix(subset_transformed_counts),
        type = "heatmap",
        colorscale = "Viridis"
      )
      hm_plt
    })
    
    # Dropdown for gene selection
    output$GeneSelect <- renderUI({
      selectInput('SelectGene', 'Select Gene of Interest', choices = rownames(dds), selected = rownames(dds)[1])
    })
    
    # Gene count plot
    output$count_plt <- renderPlotly({
      d <- plotCounts(dds, gene = input$SelectGene, intgroup = "condition", returnData = TRUE)
      count_plt <- ggplot(d, aes(x = condition, y = count, color = condition)) +
        geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
        scale_y_log10() +
        theme_bw() +
        labs(x = "Condition", y = "Count (log scale)")
      ggplotly(count_plt)
    })
  })
  
  # Placeholder text output (can be removed if unused)
  output$value <- renderText({ input$caption })
}

# Run Shiny app
shinyApp(ui = ui, server = server)
