#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(plotly)


#load data outside the server function:
#source('../code/read_data.R') #defined in global.R #https://shiny.rstudio.com/articles/scoping.html

#example app /wehisan/bioinf/Bioinformatics/SNPchipdata/haloom/brain-cox2/brain-cox2


shinyServer(function(input, output) {

### tables ###
    
    sID <- debounce(reactive({str_remove_all(input$ed_siteID," ")}),1000)
    geneID <- debounce(reactive({str_remove_all(input$ed_geneID," ")}),1000) #slow down UI reactivity
    
    ed_my_ensg     <- reactive({ ifelse(str_detect(geneID(), 'ENSG0'), geneID(),
                                     anno %>% filter(gene_name == geneID()) %>% pull(gene_id) ) })
    ed_my_symbol   <- reactive({ ifelse(!is.na(geneID()),
                                     anno %>% filter(gene_id == ed_my_ensg()) %>% pull(gene_name),
                                     'Error: Gene not found') })
    
    output$ED_ENSG_TITLE <- renderText(paste(ed_my_ensg(), ed_my_symbol(), sep=" "))
    
    count_tbl <- reactive(neuType_site_stats %>% filter(siteID==sID()) %>% select(-siteID)) 
    output$COUNT_TBL <- renderTable(count_tbl())
    
    desc_tbl <- reactive(dsc_tbl_forShiny %>% filter(`Site ID`==sID()) %>% 
                             gather("reference_db",'record')) 
    output$DESC_TBL <-renderTable(desc_tbl())
        
    scz_asd_tbl <- reactive(scz_asd_long %>% filter(siteID==sID()) %>% select(-siteID))
    output$SCZ_ASD_TBL <- renderTable(scz_asd_tbl())
    
    
### histogram ###    
    
    hstPlot <- reactive(meanAP_distn_binned %>% ggplot(aes(x=meanAP)) + 
                            geom_histogram(fill='dodgerblue') + 
                            geom_point( data = meanAP_distn_binned %>% filter(siteID==sID()),
                                        aes(x=upr,y=hght),col='red',size=4) +
                            xlab('Alternate (edited) allele frequency averaged across cells') +
                            ylab('Edited site count') +
                            theme(axis.text = element_text(size=14), axis.title=element_text(size=16)) 
    )
    output$edHist <- renderPlot(hstPlot())


### site Ed per Cell ###
    
    site_by_gp <- reactive(sites_cells_slim %>% filter(siteID == sID()) %>% 
                               ggplot(aes(x=Group_ID,y=altProp,group=neuType,col=neuType)) + 
                               geom_point(size=2, 
                                          position=position_jitterdodge(dodge.width=0.8,jitter.width = 0.1,jitter.height=0.01)  ) +
                               xlab('Neuronal subgroup') + ylab('Edited allele proportion') + labs(col='Neuronal subtype') +
                               theme(axis.text   = element_text(size=14), axis.title=element_text(size=16), 
                                     legend.text = element_text(size=16), legend.position='bottom') + labs(col='') 
                               #bottom_legend
                           )
    output$SITE_BY_GP <- renderPlot(site_by_gp())
    
    site_by_type_area <- reactive(sites_cells_slim %>% filter(siteID == sID()) %>% 
                                      ggplot(aes(x=area,y=altProp,group=neuType,group=neuType,col=neuType)) + 
                                      geom_point(size=2, 
                                                 position=position_jitterdodge(dodge.width=0.8,jitter.width = 0.15,jitter.height=0.01)  ) +
                                      xlab('Cortical region (Brodmann area)') + ylab('Edited allele proportion') + 
                                      labs(col='Neuronal subtype') +
                                      theme(legend.position = 'none', axis.text = element_text(size=14), axis.title=element_text(size=16))
                                      #no_legend
                                      
                                  
                        )
    output$SITE_BY_TYPE_AREA <- renderPlot(site_by_type_area())
    
        
    

### geiVtpm volcano ###    
    
    
gt_volc <- reactive(volcPlot_bg +
                    ggnewscale::new_scale_color() +
                    geom_point(data  = geiVtpm %>% filter(ENSG == ed_my_ensg()), 
                               aes(x = estimate,y = y_p), col='red',size=4) +
                    annotate(geom = "text", x=-0.04, y=100, label = ed_my_symbol(),cex=6,col="red") +
                    xlab('Estimate (association between gene expression and global editing rate)') +
                    ylab('-log10(p value)') +    
                    theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
                          legend.text = element_text(size=16), legend.title = element_text(size=16))
                    )
    output$GT_VOLC <- renderPlot(gt_volc())

    
    
### plotly ###
### sites FC over gene length ###   
    
    chrID <- reactive(hqEd_site_logFC_geneENR_REDI_RM %>% filter(gene_name == geneID()) %>% pull(Chr))
    
    p_annot <-reactive(
            hqEd_site_logFC_geneENR_REDI_RM %>% 
            #filter(gene_name==geneID()) %>% 
            filter(gene_id==ed_my_ensg()) %>% 
            arrange(Locus) %>% mutate(meanFC = mean(logFC), maxLoc = max(Locus),minLoc=min(Locus)) %>% 
            slice(n=1) %>% 
            mutate(rnge = maxLoc - minLoc) %>% 
            mutate(y_pos = meanFC + (0.75*(abs(meanFC))),
                   x_pos = Locus  + (0.1 * rnge) ) )
    
    p <- reactive(
        hqEd_site_logFC_geneENR_REDI_RM %>% 
        #filter(gene_name==geneID()) %>% 
        filter(gene_id==ed_my_ensg()) %>%
        select(-gene_name,-gene_id) %>% distinct() %>%  #avoid redundant datapoints
        ggplot(aes(x=Locus,y=logFC, col=REDI_RM)) + 
        #geom_point(aes(text=paste0(siteID,' ',-log10(PValue),' ',REDI_RM))) +
        #geom_point(aes(text=paste0(siteID,' \n logFC = ',logFC_rnd,' p = ',pV_sci))) +
        geom_point(aes(text = paste0(siteID,'\n ',genicFeature_summ,'\nlogFC = ',logFC_rnd,' p = ',pV_sci),
                       size = (-log10(PValue))), show.legend=FALSE) +
        scale_color_brewer(palette = 'Paired') + # ggtitle(geneID()) + 
        geom_hline(yintercept = 0,lty=2) +
        geom_segment(aes(x=min(Locus),xend=max(Locus),y=mean(logFC),yend=mean(logFC)), lwd=0.5, col='red') +
        #annotate(geom='text', aes(x = min(Locus), y= (mean(logFC) + 0.2*mean(logFC))), col='red', label ='mean logFC' ) +
        geom_text(data = p_annot(), aes(x = x_pos, y= y_pos), col='red', label ='mean logFC', hjust=1, cex=4, check_overlap = TRUE ) +
        xlab(paste0('Chr ',chrID(), ', position (BP; Hg38 build)')) +
        ylab('Log editing fold-change (Excitatory vs Inhibitory neurons) ') +
        #bottom_legend 
        theme(legend.position='bottom')
        
        
    )
    
   
     output$PLOTLY_GENE <- renderPlotly({
         validate(need( 
                nrow(hqEd_site_logFC_geneENR_REDI_RM %>% filter(gene_id ==ed_my_ensg())) >0, 
               'Sorry no variable editing sites were detected in this gene' ))
            
         print(
                ggplotly( p() , tooltip="text")) })
     
    #reactive values: copy plotly siteIDs to clip board / or to sID() [?] https://community.rstudio.com/t/dynamic-modules-and-plotly-event-data-not-working/55138/2


     
     
############# ############# ############# ############# ############# ############# ############# #############
     #############  #############   GENE EXPRESSION DATA #############  #############   #############  
############# ############# ############# ############# ############# ############# ############# #############     

    expn_geneID <- debounce(
        reactive({str_remove_all(input$input_geneID," ")}),
        1000)      #slow down reactivity https://shiny.rstudio.com/reference/shiny/1.0.4/debounce.html
    
    my_ensg     <- reactive({ ifelse(str_detect(expn_geneID(), 'ENSG0'), expn_geneID(),
                          anno %>% filter(gene_name == expn_geneID()) %>% pull(gene_id) ) })
    my_symbol   <- reactive({ ifelse(!is.na(expn_geneID()),
                                     anno %>% filter(gene_id == my_ensg()) %>% pull(gene_name),
                                     'Error: Gene not found') })

    Lake_vals <- reactive({  subset(Lake_tpm, ENSG == my_ensg() ) })

    GSE67835_vals <- reactive({GSE67835_tpm %>% dplyr::filter(ENSG == my_ensg())})

    BS_vals <- reactive({ bs_rowData %>% filter(ensembl_gene_id == my_ensg()) })

    output$ensg_title <- renderText(paste(my_ensg(), my_symbol(), sep=" "))

    
    Lake_subset <- reactive({
  			 id <- showNotification("Reading data...", duration = NULL, closeButton = FALSE)
       			 on.exit(removeNotification(id), add = TRUE)
        
        Lake_meta %>%
            left_join(Lake_vals(), by = c('SRA' = 'key')) %>%
            arrange(TPM) %>% dplyr::filter(neuType!="NoN") %>%
            mutate(cortical_area = factor(area, levels=c('BA8','BA10','BA21','BA22','BA41','BA17'))) %>%
            mutate(Neuronal_type = factor(Neuronal_type,levels=c('Excitatory','Inhibitory')))
           })
            
    
    Lake_scatter <- reactive( 
        Lake_subset() %>% 
            ggplot() +
            geom_point(data = . %>% filter(  is.na(TPM)), aes(x=dim1, y=dim2, col=log10(TPM)), cex=0.75) +
            geom_point(data = . %>% filter(! is.na(TPM)), aes(x=dim1, y=dim2, col=log10(TPM)), cex=0.75) +
            scale_color_viridis_c() +
            geom_text(data=tribble(~nType , ~X, ~Y, 'Excitatory',2,3,'Inhibitory',-10,3),
                      aes(x=X,y=Y,label=nType), fontface='bold') +
            ggtitle('' , subtitle = 'scRNA-seq of 3127 neurons from one donor') +
            theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
              legend.text = element_text(size=16), legend.title = element_text(size=16))
    )
    
    output$LAKE_SCATTER <- renderPlot( Lake_scatter() )
    

    Lake_violin_GroupID <- reactive(
        
        Lake_subset() %>% 
            ggplot(aes(x=Group_ID,y=log10(TPM))) +
            geom_violin(aes(fill=Group_ID),alpha=0.5) +
            geom_jitter(col='black',cex=0.5,height=0,width=0.1) +
            scale_fill_manual(values = group_cols) +
            facet_wrap(~Neuronal_type, scales='free_x', nrow = 1 ) +
            #no_legend +
            ggtitle('', subtitle = 'Neuronal sub-group') +
            theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
                 legend.text = element_text(size=16), legend.title = element_text(size=16), legend.position = 'none')
    )
    
    output$LAKE_VIOLIN_GPID <- renderPlot( Lake_violin_GroupID() )  

    
    Lake_violin_area <- reactive(

        Lake_subset() %>% 
            ggplot(aes(x=area,y=log10(TPM),col=Neuronal_type)) + 
            geom_jitter(position=position_jitterdodge(jitter.height = 0,jitter.width = 0.1,dodge.width = 0.8),cex=0.75) + 
            geom_boxplot(alpha=0) +
            xlab('') +
        ggtitle('', subtitle = 'Cortical region (Brodmann area)') +
            theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
                  legend.text = element_text(size=16), legend.title = element_text(size=16),
                  legend.position='bottom') +
            labs(col='')
        )

    output$LAKE_VIOLIN_AREA <- renderPlot( Lake_violin_area() )

    
    GSE67835_scatter <- reactive(

        GSE67835_plot_data %>% left_join(GSE67835_vals(), by=c('sample'='key')) %>%
            filter(!is.na(cell_type)) %>%
            # mutate(cell_type=factor(cell_type, levels=c('astrocytes','endothelial','hybrid','microglia','neurons',
            #                                             'oligodendrocytes','OPC','fetal_replicating','fetal_quiescent'))) %>%
            ggplot(aes(x=dim1,y=dim2)) +
            geom_point( alpha = 0.5, aes(col = log10(value))) +
            ggrepel::geom_label_repel(data = GSE67835_data_labels,
                                      aes(x=dim1,y=dim2, fill=cell_type, label = cell_type),
                                      nudge_y = 0.5,nudge_x = -2,
                                      segment.alpha = 0, alpha=0.25, show.legend = FALSE, seed = 1234) +
            ggrepel::geom_label_repel(data = GSE67835_data_labels,
                                      aes(x=dim1,y=dim2, fill=cell_type, label = cell_type),
                                      nudge_y = 0.5,nudge_x = -2,
                                      label.size=NA, segment.alpha=0, alpha=1, fill=NA, show.legend = FALSE, seed=1234) +
            scale_color_viridis_c() +
            scale_fill_manual(values = GSE67835_cols) +
            ggtitle('', subtitle = '365 brain cells from 12 donors') +
            theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
                  legend.text = element_text(size=16), legend.title = element_text(size=16))
        
    )
    
    output$GSE67835_SCATTER <- renderPlot(  GSE67835_scatter() )
    

    GSE67835_violin <- reactive(

        GSE67835_plot_data %>% left_join( GSE67835_vals(),by=c('sample'='key')) %>%
            filter(!is.na(cell_type)) %>%
             mutate(cell_type=factor(cell_type, levels=c('astrocytes','endothelial','hybrid','microglia','neurons',
                                                         'oligodendrocytes','OPC','fetal_replicating','fetal_quiescent'))) %>%
            ggplot(aes(x=cell_type, y=log10(value))) +
            geom_violin(aes(fill=cell_type), alpha=0.25, outlier.alpha = 0, show.legend = FALSE) +
            geom_jitter(width=0.15, height=0, cex=0.75, show.legend = FALSE) +
            scale_fill_manual(values = GSE67835_cols) +
            ylab('log10(transcripts per million)') + xlab("") +
            theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
                  legend.text = element_text(size=16), legend.title = element_text(size=16))
    )
    
    output$GSE67835_VIOLIN <- renderPlot( GSE67835_violin() )

    

#     
#     output$GSE84465_scatter <- renderPlot({
#         
#         GSE84465_vals_r <- GSE84465_vals()
#         
#         GSE84465_plot_data %>% left_join(GSE84465_vals_r, by=c('sample'='key')) %>% 
#             mutate(status = case_when(facet_groups == 'Immune cell' & dim1 < 4 ~ 'drop',
#                                          facet_groups=='OPC' & dim1< (-11) ~ 'drop',
#                                          facet_groups=='Neoplastic' & dim1 > -5  ~ 'drop',
#                                          facet_groups=="Neoplastic" & dim1 < (-11) ~ 'drop',
#                                          facet_groups=='other'& dim1>0 ~ 'drop',
#                                          TRUE ~ 'retain')) %>% 
#             filter(status=='retain') %>% 
#             mutate(facet_groups = case_when(facet_groups =='Neoplastic' & dim2 < 5 ~ 'Neoplastic_1',
#                                             facet_groups =='Neoplastic' & dim2 > 5 ~ 'Neoplastic_2',
#                                             TRUE ~ facet_groups)) %>% 
#             #mutate(clust_split = ifelse(dim1 < 5, 1, 2)) %>% 
#             ggplot(aes(x=dim1,y=dim2)) + 
#             geom_point( alpha = 0.5, cex = 1.2, aes(col = log10(value))) +
#             ggrepel::geom_label_repel(data = GSE84465_data_labels %>% filter(!str_detect( cell_type,'Neoplastic|Immune|OPC')) , #%>% mutate(clust_split=ifelse(dim1<5,1,2)), 
#                                       aes(x=dim1,y=dim2, fill=cell_type, label = cell_type), 
#                                       nudge_y = 1.5, segment.alpha = 0, alpha=0.25, show.legend = FALSE, seed = 1234) +
#             ggrepel::geom_label_repel(data = GSE84465_data_labels %>% filter(!str_detect( cell_type,'Neoplastic|Immune|OPC')), #%>% mutate(clust_split=ifelse(dim1<5,1,2)), 
#                                       aes(x=dim1,y=dim2, label = cell_type), 
#                                       nudge_y = 1.5, label.size=NA, segment.alpha=0, alpha=1, fill=NA, show.legend = FALSE, seed=1234) +
#             scale_color_viridis_c() +
#             scale_fill_manual(values = GSE84465_cols) +
#             facet_wrap(~ facet_groups, scales='free', ncol=2, dir = 'v') +
#             ggtitle('Darmanis et al Cell Reports 2017', subtitle = '2852 cells from 4 glioblastoma patients')
#         
#             })
#     
#     output$GSE84465_violin <- renderPlot({
#         
#        GSE84465_vals_r <- GSE84465_vals()
#        
#        GSE84465_plot_data %>% left_join(GSE84465_vals_r,by=c('sample'='key')) %>%
#             ggplot(aes(x=cell_type,y=log10(value))) + 
#             geom_violin(aes(fill=cell_type), alpha=0.25, outlier.alpha = 0, show.legend = FALSE) + 
#             geom_jitter(width=0.15, height=0, cex=0.75, show.legend = FALSE) +
#             scale_fill_manual(values = GSE84465_cols) +
#             ylab('log10(transcripts per million)') + xlab("") 
#     })
#     

    # BrainSpan_box <- reactive(
    # 
    #     gene_brainSp_data <- left_join(BS_vals(), bs_expn, by=c('row_num'='X1')) %>%
    #         #filter(gene_symbol == mygene ) %>%
    #         gather(key,value,starts_with('X')) %>% mutate(key=as.numeric(str_remove(key,"X"))) %>%
    #         left_join(bs_colData, by=c("key" = "column_num")) %>%
    #         mutate(age_bin=cut_number(age,5)) %>% filter(!is.na(age_bin))
    # 
    #     gene_brainSp_data %>%
    #         ggplot(aes(x=age_bin, col=gene_symbol, fill=gene_symbol, y=log10(value))) +
    #         geom_boxplot(aes(), position='dodge',alpha=0.25, outlier.colour = NULL) +
    #         geom_point(position=position_jitterdodge(dodge.width = 0.8,jitter.width = 0.1, jitter.height = 0), cex=0.25, col='black') +
    #         stat_summary(fun.y=median, geom="line",  aes(group=gene_symbol), position = position_dodge(width=0.8)) +
    #         stat_summary(fun.y=median, geom="point", aes(group=gene_symbol), position = position_dodge(width=0.8)) +
    #         facet_wrap(~ lobe, scales='free_y', ncol = 2, dir = "v") +
    #         xlab('Age group (years)') + ylab('log10(length-normalized transcription)') +
    #         theme(legend.position='bottom') +
    #         ggtitle('BrainSpan', subtitle = 'Bulk RNA-seq from 42 individuals grouped by brain lobe') +
    #         theme(axis.text = element_text(size=14), axis.title=element_text(size=16),
    #               legend.text = element_text(size=16), legend.title = element_text(size=16))
    # )
    # output$BRAINSPAN_BOX <- renderPlot(BrainSpan_box())


    
    outputOptions(output, "LAKE_SCATTER", suspendWhenHidden = FALSE)
    outputOptions(output, "LAKE_VIOLIN_GPID", suspendWhenHidden = FALSE)
    outputOptions(output, "LAKE_VIOLIN_AREA", suspendWhenHidden = FALSE)
    outputOptions(output, "GSE67835_SCATTER", suspendWhenHidden = FALSE)
    outputOptions(output, "GSE67835_VIOLIN", suspendWhenHidden = FALSE)
    #outputOptions(output, "BRAINSPAN_BOX", suspendWhenHidden = FALSE)
    #outputOptions(output, "GSE84465_scatter", suspendWhenHidden = FALSE)
    #outputOptions(output, "GSE84465_violin", suspendWhenHidden = FALSE)
    #outputOptions(output, "", suspendWhenHidden = FALSE)


})






