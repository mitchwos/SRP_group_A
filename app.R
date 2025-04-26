#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

#uploading relevant library 
library(shiny)
library(bslib)

# ---- UI ----
ui <- navbarPage(
  title = "Darmanis Reanalysis",
  id = "mainNav",
  theme = bs_theme(version = 5, bootswatch = "darkly"), #Theme of the website using bootstrap
  
  #Home Tab
  tabPanel("Home",
           h1("Revisiting the Darmanis et al., Single-Cell Brain Atlas", class = "text-primary"),
           h4("Welcome to our website where we offer an in-depth re-analysis of the single-cell brain atlas produced by Darmanis et al."), #Link to the study 
           h3("Background Information:"),
           p("This study was conducted to gain a deeper understanding of the molecular diversity of human brain cells at a single-cell level. Researchers used single-cell RNA sequencing to analyse 466 individual cells from both adult and fetal human brains. Darmanis et al. sought to explore and classify major neuronal, glial, and vascular cell types found within the human brain."),
           p("For an in-depth explanation of this study, check out: ",
             a("Darmanis et al. Single-Cell Brain Atlas study", href = "https://pubmed.ncbi.nlm.nih.gov/26060301/", target = "_blank"))
  ),
  
  #Original Pipeline Tab 
  tabPanel("Original Pipeline",
           h1("Recreating the Darmanis et al. Pipeline", class = "text-primary"),
           h3("Methods used"),
           tags$ul(
             tags$li("Gene count and metadata from ", a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835", target = "_blank", "GSE67835")), #link to GEO
             tags$li("Analysis in R (2025) using default parameters"),
             tags$li("Pairwise distances via `scde`, differential expression also assumed from `scde`"),
             tags$li("Metadata retrieved using GEOquery and integrated using `SingleCellExperiment`"),
             tags$li("Dimensionality reduction with `Rtsne`"),
             tags$li("Unbiased clustering using `mclust` with BIC + GMMs and EM algorithm"),
             tags$li("Markers identified with PanglaoDB"),
             tags$li("Biased clustering based on gene selection from Darmanis et al. Fig 1B"),
             tags$li("Visualisation using `scatterplot3d`, `pheatmap`, and `ggplot2`")
           ),
           h3("Challenges faced"),
           p("There were several difficulties during the reanalysis. The original study didn’t give enough detail about the exact methods, software, or settings they used, which made it hard to follow their approach closely. Some data, like the mouse brain RNA-seq reference, was no longer available, so we had to guess the genes used for certain parts of the analysis by looking at their figures. We also found small differences in clustering, such as getting three fetal clusters instead of two, and some clusters looked very similar to each other. Some tools had to be downgraded to work properly, which took extra time. Because of all this, certain parts of the original study—like detailed comparisons between fetal and adult neurons—couldn’t be fully recreated. However, most of our results still matched the original findings quite well.
"),
           h3("Results"),
           p("Clustering of 466 brain cells revealed 10 distinct groups, with adult cells in clusters 1–7 and fetal cells in clusters 8–10. Using Gaussian Mixture Models and gene markers, adult clusters aligned with cell types like neurons, astrocytes, microglia, and oligodendrocytes, while fetal clusters lacked definitive markers. Biased clustering refined these results and showed agreement with the unbiased method, although some neuron clusters overlapped. Analysis of MHC-I genes indicated immune activity in adult neurons, whereas fetal cells showed limited expression, primarily in replicating subtypes. These findings highlight distinct developmental roles, with fetal cells emphasizing neurogenesis and adult cells specializing in immune functions."),
           div(
             style ="display: flex; flex-direction: column; align-items: flex-start;", # allows the picture to stay on the left of the webpage
             tags$p(class = "fw-bold", "Click the image to enlarge"),
             tags$a(href = "BIC.png", target = "_blank", tags$img(src = "BIC.png", width = "30%")),
             tags$p("10 components by BIC for GMMs. Unbiased brain cell clusters."),
             tags$a(href = "3D.png", target = "_blank", tags$img(src = "3D.png", width = "30%")),
             tags$p("3D clustering of single cells."),
             tags$a(href = "UTF-8heatmap.png", target = "_blank", tags$img(src = "UTF-8heatmap.png", width = "20%")),
             tags$p("Heatmap of gene expression by cluster."),
             tags$p("Orange: oligodendrocytes, Yellow: astrocytes, Green: microglia, Pink: endothelial, Blue: OPC, Red: neurons"),
             tags$a(href = "cell_type_assignment.png", target = "_blank", tags$img(src = "cell_type_assignment.png", width = "30%")),
             tags$p("Comparison of cluster-based vs known assignments."),
             tags$a(href = "MHC_I.png", target = "_blank", tags$img(src = "MHC_I.png", width = "30%")),
             tags$p("MHCI pathway expression in adult and fetal cells.")
           )
  ),
  
  #Alternative Pipeline Tab 
  tabPanel("Alternative Pipeline",
           h1("Alternative Pipeline", class = "text-primary"),
           h3("Methods"),
           tags$ul(
             tags$li("Gene count matrix from HISAT2 + featureCounts"),
             tags$li("Normalization with LogNormalize (scale 1e6), HVG selection"),
             tags$li("Dimensionality reduction using PCA + UMAP, clustering via Louvain"),
             tags$li("Gene symbols mapped using biomaRt"),
             tags$li("Metadata from GEO:", a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835", "GSE67835")), #link to GEO
             tags$li("Markers found with `FindAllMarkers`, differential expression across conditions"),
             tags$li("Plots made with ggplot2, dplyr, patchwork")
           ),
           h3("Challenges faced"),
           p("A key challenge faced during the analysis was the need to convert Ensembl gene identifiers to HGNC gene symbols, as the Ensembl IDs were difficult to match with the metadata and hindered interpretation of results. Additionally, the lack of detailed methods provided by the original authors made it challenging to replicate results, especially since the use of an alternative pipeline involved different methodologies. Therefore, the absence of clear guidance further complicated the analysis process."),
           h3("Results"),
           p("The alternative pipeline identified 10 clusters using the UMAP projection, aligning cell groups with distinct biological identites. The adults clusters aligned to cell types like neurons, astrocytes, microglia and oligondendrocytes, while fetal clusters reflected development stages including quiescent and replicating states.Marker gene analysis highlighted unique roles such as immune responses and neurodevelopment. Differential gene expression showed fetal cells focusing on neurogenesis and proliferation, while adult cells exhibited specialized functions like myelination and signaling. MHC-I gene expression was minimal in fetal cells but strongly expressed in adult immune-related cells"),
           div(
             style ="display: flex; flex-direction: column; align-items: flex-start;", # allows the picture to stay on the left of the webpage
             tags$p(class = "fw-bold", "Click the image to enlarge"),
             tags$a(href = "umap.png", target = "_blank", tags$img(src = "umap.png", width = "25%")),
             tags$p("UMAP visualisation of clustering and developmental stages.  (A) Unbiased clustering identifies ten distinct clusters labeled 1–10. (B) Developmental stages of cells categorized as Adult, Fetal_Quiescent, and Fetal_Replicating. (C) Biased clustering assigns cell types: astrocytes, endothelial cells, fetal_quiescent cells, fetal_replicating cells, hybrid cells, microglia, neurons, oligodendrocytes, and OPCs."),
             tags$a(href = "unbiased.png", target = "_blank", tags$img(src = "unbiased.png", width = "25%")),
             tags$p("Marker gene expression across unbiased clusters.  The heatmap illustrates the top five marker genes for each of the 10 unbiased clusters. Columns represent clusters, while rows denote genes. Colour intensity reflects expression levels, with yellow indicating high expression, black medium, and purple low. Cluster identities are colour-coded above the heatmap"),
             tags$a(href = "biased.png", target = "_blank", tags$img(src = "biased.png", width = "25%")),
             tags$p("Heatmap of marker gene expression across biased clusters.  The heatmap displays the top five marker genes for each biased cell type, including astrocytes, endothelial cells, fetal_quiescent cells, fetal_replicating cells, hybrid cells, microglia, neurons, oligodendrocytes, and OPCs. Rows represent marker genes, while columns denote cell types. Gene expression levels are colour-coded, with yellow indicating high expression and black indicating low expression."),
             tags$a(href = "differential.png", target = "_blank", tags$img(src = "differential.png", width = "25%")),
             tags$p("Differential gene expression between fetal and adult cell types.  The heatmap depicts gene expression levels of the top differentially expressed genes across various cell types, including astrocytes, endothelial cells, hybrid cells, microglia, neurons, oligodendrocytes, OPCs, fetal_quiescent cells, and fetal_replicating cells. Yellow represents high expression, black medium, and purple low. Gene names are listed on the left, while cell types are color-coded at the top."),
             tags$a(href = "MHC1.png", target = "_blank", tags$img(src = "MHC1.png", width = "25%")),
           tags$p("MHC1 gene expression by cell type.  The dot plot illustrates the expression of MHC1 pathway genes (HLA-A, HLA-B, HLA-C, TAPBP, CALR, ERAP1, B2M, PDI4A, HSPAS) across various cell types, including OPC, oligodendrocytes, neurons, microglia, hybrid cells, fetal_quiescent, fetal_replicating, endothelial, and astrocytes. Dot colour represents the average gene expression (blue: high expression; white: low expression; purple: negative expression), while dot size reflects the percentage of cells expressing the gene.")
             )
  ),
  
  # Cluster Exploration Tab 
  tabPanel("Cluster Exploration",
           sidebarLayout(
             sidebarPanel(
               # Radio buttons to click biased and unbiased
               radioButtons("search_type", "Choose Search Type:", choices = c("Unbiased Clusters" = "unbiased", "Biased Clusters" = "biased")),
               # Drop down for unbiased 
               conditionalPanel(
                 condition = "input.search_type == 'unbiased'",
                 selectInput("unbiased_cluster", "Select Cluster (1-10):",
                             choices = paste("Cluster", 1:10))
               ),
               # Shows the drop down for biased 
               conditionalPanel(
                 condition = "input.search_type == 'biased'",
                 selectInput("biased_type", "Select Cell Type:",
                             choices = c("Astrocytes", "Neurons", "Microglia", "Oligodendrocytes", "OPCs", "Endothelial Cells", 
                                         "Fetal replicating", "Hybrid", "Fetal quiescent"))
               )
             ),
             mainPanel(
               h4("Results"), #Heading for results
               uiOutput("result"), #Place holder where the genes will be displayed
               uiOutput("image")  #place holder for images 
             )
           )
  )
)

#  SERVER 
server <- function(input, output, session) {
  unbiased_genes <- list(                           #Lists generated to add into the interactive part
    "Cluster 1" = c("PLAG2GS", "CHI3L1", "MT1H"),
    "Cluster 2" = c("KLHL1", "ZNF93", "DRAXIN", "FEZF2", "NKAIN1"),
    "Cluster 3" = c("NDNF", "EREG", "PTGFR", "PRDX1P1", "EGFL6", "CACNA1G", "IQSEC3-AS1", "COL19A1"),
    "Cluster 4" = c("SOX11", "MR041", "DRAXIN"),
    "Cluster 5" = c("EGFL6", "CACNA1G", "IQSEC3-AS1", "STRC"),
    "Cluster 6" = c("TBX3", "IL6", "CRISPLD2", "LIMS2"),
    "Cluster 7" = c("GPR17", "LIMS2", "CLCA4", "KLK6", "ACSL3P1", "LINC02882"),
    "Cluster 8" = c("CLEC4G", "GPR88", "NRGN", "C3orf80"),
    "Cluster 9" = c("C1QA", "LPAR5", "FOLR2", "SIGLEC8"),
    "Cluster 10" = c("RRM2", "E2F2", "E2F8", "CDK1", "HJURP")
  )
  
  biased_genes <- list(
    "Astrocytes" = c("HN2F2-AS1", "RNU6-36SP", "CH13L1", "LIN00391", "SLC39A12-AS1"),    #Lists generated to add into the interactive part
    "Neurons" = c("ARHGD1G", "FNBCQ", "KCNJ9", "COR06"),
    "Microglia" = c("CD300E", "LPAR5", "LINC01480", "HPGP5"),
    "Endothelial Cells" = c("IGF2", "TWIST1", "TBX18", "SEMA3G"),
    "Fetal replicating" = c("AURKB", "RRM2", "TCP2A", "HJURP", "SKA1"),
    "Hybrid" = c("KH5RP", "MAGEB17", "A1DH7", "NR2F2-AS1"),
    "OPCs" = c("DISC1-IT1", "MROH9", "STON1"),
    "Oligodendrocytes" = c("ACSL3P1", "CLA4", "KLK6", "SLC5A11", "MAG"),
    "Fetal quiescent" = c("DRAXIN", "ADRA2A", "SCUB1", "NEUROD6")
  )
  
  output$result <- renderUI({               #Allows the output depending on user selection 
    if (input$search_type == "unbiased") {
      genes <- unbiased_genes[[input$unbiased_cluster]]
      tagList(
        h4(paste("Unbiased Cluster:", input$unbiased_cluster)),
        tags$ul(lapply(genes, tags$li))
      )
    } else {
      genes <- biased_genes[[input$biased_type]]
      tagList(
        h4(paste("Biased Cluster: Cell Type -", input$biased_type)),
        tags$ul(lapply(genes, tags$li))
      )
    }
  })
  
  output$image <- renderUI({              #Pictures added and depending on user selection 
    if (input$search_type == "unbiased") {
      tags$img(src = "unbiased.png", alt = "Unbiased Clusters", style = "max-width:100%; border:1px solid black;")
    } else {
      tags$img(src = "biased.png", alt = "Biased Clusters", style = "max-width:100%; border:1px solid black;")
    }
  })
}

# ---- Run App ----
shinyApp(ui = ui, server = server)

