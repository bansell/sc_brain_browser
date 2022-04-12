#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)


# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Human CNS gene expression & RNA editing browser"),

            tabsetPanel(type = 'tabs', #https://shiny.rstudio.com/articles/tabsets.html
            tabPanel('About',
                     br(),
                     h5("Welcome to the CNS browser. Created by Brendan Ansell, Bahlo Laboratory, Walter + Eliza Hall Institute"),
                     br(),
                     # h5("Please use the Browse tab to search for your gene of interest"),
                     # br(),
                     h4('Abstracts of represented publications'),
                     br(),
                     h5(a("Lake et al Science 2016", href='https://science.sciencemag.org/content/352/6293/1586'),
                        br(),br(),
                        "The human brain has enormously complex cellular diversity and connectivities fundamental to our neural functions, yet difficulties in interrogating individual neurons has impeded understanding of the underlying transcriptional landscape. We developed a scalable approach to sequence and quantify RNA molecules in isolated neuronal nuclei from a postmortem brain, generating 3227 sets of single-neuron data from six distinct regions of the cerebral cortex. Using an iterative clustering and classification approach, we identified 16 neuronal subtypes that were further annotated on the basis of known markers and cortical cytoarchitecture. These data demonstrate a robust and scalable method for identifying and categorizing single nuclear transcriptomes, revealing shared genes sufficient to distinguish previously unknown and orthologous neuronal subtypes as well as regional identity and transcriptomic heterogeneity within the human brain."),
                     br(),
                     img(src="NeuType_forShiny.png",height=500,width=800),
                     br(),
                     h5(a("Darmanis et al PNAS 2015", href="https://www.ncbi.nlm.nih.gov/pubmed/26060301"),
                        br(),br(),
                        "The human brain is a tissue of vast complexity in terms of the cell types it comprises. Conventional approaches to classifying cell types in the human brain at single cell resolution have been limited to exploring relatively few markers and therefore have provided a limited molecular characterization of any given cell type. We used single cell RNA sequencing on 466 cells to capture the cellular complexity of the adult and fetal human brain at a whole transcriptome level. Healthy adult temporal lobe tissue was obtained during surgical procedures where otherwise normal tissue was removed to gain access to deeper hippocampal pathology in patients with medical refractory seizures. We were able to classify individual cells into all of the major neuronal, glial, and vascular cell types in the brain. We were able to divide neurons into individual communities and show that these communities preserve the categorization of interneuron subtypes that is typically observed with the use of classic interneuron markers. We then used single cell RNA sequencing on fetal human cortical neurons to identify genes that are differentially expressed between fetal and adult neurons and those genes that display an expression gradient that reflects the transition between replicating and quiescent fetal neuronal populations. Finally, we observed the expression of major histocompatibility complex type I genes in a subset of adult neurons, but not fetal neurons. The work presented here demonstrates the applicability of single cell RNA sequencing on the study of the adult human brain and constitutes a first step toward a comprehensive cellular atlas of the human brain."),
                     br(),
                     h5(a("Darmanis et al Cell Reports 2017", href="https://www.ncbi.nlm.nih.gov/pubmed/29091775"),
                        br(),br(),
                        "Glioblastoma (GBM) is the most common primary brain cancer in adults and is notoriously difficult to treat because of its diffuse nature. We performed single-cell RNA sequencing (RNA-seq) on 3,589 cells in a cohort of four patients. We obtained cells from the tumor core as well as surrounding peripheral tissue. Our analysis revealed cellular variation in the tumor's genome and transcriptome. We were also able to identify infiltrating neoplastic cells in regions peripheral to the core lesions. Despite the existence of significant heterogeneity among neoplastic cells, we found that infiltrating GBM cells share a consistent gene signature between patients, suggesting a common mechanism of infiltration. Additionally, in investigating the immunological response to the tumors, we found transcriptionally distinct myeloid cell populations residing in the tumor core and the surrounding peritumoral space. Our data provide a detailed dissection of GBM cell types, revealing an abundance of information about tumor formation and migration."),
                     br(),
                #     h5(a("BrainSpan Database", href="http://help.brain-map.org/download/attachments/3506181/Transcriptome_Profiling.pdf?version=1&modificationDate=1382036562736&api=v2"),
                        br(),br() #,
                #        "BrainSpan, an atlas of the developing human brain, is designed as a foundational resource for studying transcriptional mechanisms involved in human brain development. It is the outcome of three ARRA-funded grants through the National Institutes of Health supporting a consortium consisting of the Allen Institute for Brain Science; Yale University (Nenad Sestan, Mark B. Gerstein); the Zilkha Neurogenetic Institute of the Keck School of Medicine of the University of Southern California (James A. Knowles, Pat Levitt); the Athinoula A. Martinos Center at Massachusetts General Hospital/Harvard Medical School and MIT HST/CSAIL (Bruce Fischl); the University of California, Los Angeles (Daniel H. Geschwind); and the University of Texas Southwestern Medical Center (Hao Huang) with strong collaborative support from the Genes, Cognition and Psychosis Program, which is part of the Intramural Research Program of NIMH, NIH (Thomas M. Hyde, Joel E. Kleinman, Daniel R. Weinberger). All data are publicly accessible via the ALLEN BRAIN ATLAS data portal or directly at www.developinghumanbrain.org.")
            ),
           
            tabPanel('Single cell gene expression',
            sidebarLayout(
                sidebarPanel(textInput('input_geneID', 'ENSEMBL Gene ID or Gene Symbol:', value = 'ENSG00000160710'),
                             h5('Please allow a few seconds for your data to load')),
                mainPanel(
             h3(textOutput('ensg_title')) 
             ,br()
             ,h4('Healthy brain single cell gene expression')
             ,br()
             ,h4('Lake et al, Science 2016')
            ,plotOutput('LAKE_SCATTER',width=1000)
            ,plotOutput('LAKE_VIOLIN_GPID',width=1000)
            ,plotOutput('LAKE_VIOLIN_AREA',width=1000)
            ,br(),br()
            ,h4('Darmanis et al, PNAS 2015')
            ,plotOutput('GSE67835_SCATTER',width=1100)
            ,plotOutput('GSE67835_VIOLIN',width=1000)
            #,plotOutput('BrainSpan_box', height="900px", width="1000px")
                ))),

            
            tabPanel('RNA editing',
                     #sidebarLayout(
                     #    sidebarPanel(
                      mainPanel(
                          fluidRow(
                            br()
                            ,h4("Welcome to the CNS browser for healthy brain single-neuron RNA editing. Here you can explore RNA editing sites quantified in transcriptomes from single neuronal nuclei published by", a('Lake et al, Science 2016', href='https://science.sciencemag.org/content/352/6293/1586'),", integrated with multiple RNA editing datasets.") 
                            # and ", a('Darmanis et al, PNAS 2015', href='https://www.ncbi.nlm.nih.gov/pubmed/26060301'),".")
                            ,br()
                            ,h5('References for the datasets used in this server are provided at the bottom of the page.')
                            ,br()
                            ,h5('Enter your editing site of interest (Hg38 coordinates) in the box below, or scroll down to investigate editing across entire genes.')
                            ,br()
                            ,textInput('ed_siteID', 'RNA editing site ID (chr_position):', value = '1_39682184')
                              ,tableOutput('COUNT_TBL')
                              ,splitLayout(tableOutput('DESC_TBL'), tableOutput('SCZ_ASD_TBL'))
                              ,plotOutput('edHist', width = 900)
                              ,br(),br()
                              ,h5('Edited allele proportion per cell')
                              ,h6('Grouped by neuronal subgroup')
                              ,plotOutput('SITE_BY_GP',   width = 900)   
                              ,br(),br()
                              ,h6('Grouped by cortical region')
                              ,plotOutput('SITE_BY_TYPE_AREA',width = 900)
                              ,br(),br()
                              ,textInput('ed_geneID', 'Edited ENSEMBL Gene ID or Gene Symbol:', value = 'KCNIP4')  
                              ,h4(textOutput('ED_ENSG_TITLE'))
                              #,plotOutput('gene_fc_plot', width = 900)
                              ,plotlyOutput('PLOTLY_GENE', width=1200,height=600)
                              ,h5('Scroll over points to reveal editing site details. Click and drag to zoom in to a region of the plot.')
                              ,br(),br(),br()
                             ,plotOutput('GT_VOLC',     width = 900) 
                            ,br()
                            ,br()
                            ,h4('Reference datasets and publications')
                            ,a('Darmanis et al, PNAS 2015', href='https://www.ncbi.nlm.nih.gov/pubmed/26060301')
                            ,br()
                            ,a('Breen et al Nature Neuroscience 2019', href='http://dx.doi.org/10.1038/s41593-019-0463-7')
                            ,br()
                            ,a('Tran et al Nature Neuroscience 2019', href='http://dx.doi.org/10.1038/s41593-018-0287-x')
                            ,br()
                            ,a('Picardi et al Nucleic Acids Research 2017',href='http://dx.doi.org/10.1093/nar/gkw767')
                            ,br()
                            ,a('Landrum et al Nucleic Acids Research 2018',href='http://dx.doi.org/10.1093/nar/gkx1153')
                            ,br(),br(),br()
                            
                            
                                    )))
                         
            
            
            # , tabPanel('Glioblastoma',
            #          sidebarLayout(
            #            sidebarPanel(textInput('input_geneID', 'ENSEMBL Gene ID or Gene Symbol:', value = 'ENSG00000160710')),
            #            mainPanel())),
            ## sidebarLayout as above
            # ,h3(textOutput('ensg_title') )
            # ,br()
            # ,h4('Glioblastoma')
            # ,plotOutput('GSE84465_scatter', height="900px")
            # ,plotOutput('GSE84465_violin'))
                         
                         
        ))
    )



