#import data


library(here)
detach("package:here", unload=TRUE) # https://github.com/r-lib/here/issues/28#issuecomment-574271127

setwd("/home/ansell.b/sc_brain_browser_data")

library(here)
library(shiny)
library(tidyverse)
library(plotly)

#for shiny.wehi.edu.au server:

#setwd("/home/ansell.b/ShinyApps/shiny_CNS") 

#load("data/ShinyCNS.Rdata")


source('/home/ansell.b/sc_brain_browser_data/code/tidyExt.R')

source('/home/ansell.b/sc_brain_browser_data/code/colour_palettes.R')

bottom_legend <- theme(legend.position='bottom')

no_legend <- theme(legend.position = 'none')


group_cols <- c(drsimonj_pal('hot')(8),drsimonj_pal('cool')(8))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

anno <- vroom::vroom(paste0(here(),'/data/geneID_gene_name.tsv'))

#Blue Lake et al
Lake_meta <- readRDS(paste0(here(), '/data/Lake_2016/BL_metadata_highQualCells.Rds')) %>% filter(neuType!="NoN") %>%
  mutate(Neuronal_type = factor(ifelse(neuType == 'Ex','Excitatory','Inhibitory'), levels=c('Inhibitory','Excitatory')))

Lake_tpm <- vroom::vroom(paste0(here(), '/data/Lake_2016/lake_tpm_filt.tsv')) #only genes with TPM>=0.1

dsc_tbl_forShiny <- vroom::vroom(paste0(here(),'/data/Lake_2016/dt_siteStats_2020_forShiny.tsv'), delim = '\t',
                                 col_types = cols(.default = "c"))

scz_asd_long <- vroom::vroom(paste0(here(),'/data/Lake_2016/scz_asd_long.tsv'), col_types=c('ccdd')) #col_types=cols(.default='c'))
scz_asd_regions <- vroom::vroom(paste0(here(),'/data/Lake_2016/scz_asd_regions.tsv'), delim='\t')

#count_ed_tbl <- vroom::vroom(paste0(here(),'/data/Lake_2016/count_ed_tbl.tsv'), col_types=c('cciii'))
neuType_site_stats <- vroom::vroom(paste0(here(),'/data/Lake_2016/neuType_site_stats.tsv'),col_types=c('cciiiiidd'))



#for Plots
hqEd_site_logFC_geneENR_REDI_RM <- vroom::vroom(paste0(here(),'/data/Lake_2016/hqEd_site_logFC_geneENR_REDI_RM.tsv')) 

geiVtpm <- vroom::vroom(paste0(here(), '/data/Lake_2016/geiVtpm.tsv')) %>% mutate(fdr_sig=factor(fdr_sig))

volcPlot_bg <- geiVtpm %>% 
  mutate(fdr_sig= factor(ifelse(fdr_sig==0,'No','Yes'),levels=c('No','Yes'))) %>%  
  arrange(fdr_sig) %>% 
  ggplot(aes(x=estimate,y=y_p)) + 
  geom_point(aes(col=fdr_sig)) + 
  scale_color_manual(values=c('dark grey','black')) +
  xlim(-0.05,0.05) + geom_vline(xintercept = 0,lty=2) +
  labs(col='FDR < 0.05')


sites_cells_slim <- vroom::vroom(paste0(here(), '/data/Lake_2016/sites_cells_slim.tsv')) %>%
 mutate(area=factor(area, levels=paste0('BA',c(8,10,41,21,22,17)))) 


meanAP_distn_binned <-vroom::vroom(paste0(here(), '/data/Lake_2016/meanAP_distn_binned.tsv'))


########## ########## ########## ########## ########## ########## ########## ##########
#DARMANIS gene expression 
########## ########## ########## ########## ########## ########## ########## ##########

GSE67835_tpm <- vroom::vroom(paste0(here(), '/data/GSE67835_Darmanis/GSE67835_tpmSCE.tsv'))
GSE67835_plot_data <- vroom::vroom(paste0(here(), '/data/GSE67835_Darmanis/GSE67835_plot_data.tsv'))
GSE67835_data_labels <- readRDS(paste0(here(), '/data/GSE67835_Darmanis/GSE67835_data_labels.Rds'))

GSE67835_cols <- gg_color_hue(9)


########## ########## ########## ########## ########## ########## ########## ##########
#BrainSpan gene expression 
########## ########## ########## ########## ########## ########## ########## ##########

 
# bs_expn <- vroom::vroom(paste0(here(), '/data/BrainSpan/expnData.tsv'))
# bs_colData <- vroom::vroom(paste0(here(), '/data/BrainSpan/colData_cleaned.tsv'))
# bs_rowData <- vroom::vroom(paste0(here(), '/data/BrainSpan/rowData.tsv'))



