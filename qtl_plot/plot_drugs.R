library(stringr)
library(purrr)
library(rols)
library(tidyverse)
library(ggh4x)

##############
# Variables
###############

server <- 'farm'
# server <- 'openstack'

if (server == 'farm'){
  repo.dir <- '~/eqtl/code/IBDVerse-sc-eQTL-code/'
  coloc.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/2025_03_IBDverse_coloc_all_gwas/collapsed/'
  out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/therapeutics/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  coloc.dir <- paste0('/home/ubuntu/eqtl/data/IBDverse/multi_tissue/coloc-base')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/therapeutics/'
}


data.dir <- paste0(repo.dir,'/data/')

# SAVE
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE

##############
# Functions
##############

source(paste0(repo.dir,'/qtl_plot/helper_functions.R'))

##############
# Main
###############


# Create a single, reusable OlsClient
mondo <- Ontology("mondo")
efo <- Ontology("efo")
orphanet <- Ontology("ordo") 
orphanet.terms <- Terms(orphanet) # Broken for some reason so hacky fix
hp <- Ontology("HP")

otar.ibd.colocs <- read_otar(repo.dir)

# monogenic <- read_csv(paste0(data.dir, 'list_monogenic_ibd_list.csv')

up_regulators <- c("ACTIVATOR", "AGONIST")

down_regulators <- c("INHIBITOR", 
                     "ANTISENSE INHIBITOR", 
                     "ANTAGONIST", 
                     "NEGATIVE ALLOSTERIC MODULATOR", 
                     "BLOCKER")

unclear_actions <- c("MODULATOR", 
                     "STABILISER", 
                     "BINDING AGENT", 
                     "OTHER", 
                     "VACCINE ANTIGEN")

drug.data <- read_tsv(paste0(data.dir, 'coloc_chembl_moa.tsv')) %>% 
  left_join(grch38,by = c('targetId' = 'ensgene')) %>% 
  mutate(action_summary = case_when(actionType %in% up_regulators ~ 1,
                                    actionType %in% down_regulators ~ -1,
                                    TRUE ~ 0))
  # filter(action_summary != 0) %>% 

ontology.list.diseases <- lapply(unique(drug.data$diseaseId), classify_disease)
ontology.diseases.df <- as.data.frame(do.call(rbind, ontology.list.diseases)) %>% 
  mutate(V2 = str_to_sentence(V2),
         V3 = str_to_sentence(V3))
names(ontology.diseases.df) <- c('diseaseId','Disease_class','Disease_label') 


summarised.drug.data <- drug.data %>% 
  left_join(ontology.diseases.df, by='diseaseId') %>% 
  group_by(targetId, Disease_class) %>% 
  slice_max(clinicalPhase,with_ties = FALSE)

plot.summarised.drug.data <- summarised.drug.data %>% 
  mutate(phase_image = paste0(paste0(data.dir, 'semicircles/phase_',
                              ceiling(clinicalPhase),'.png')),
         is_approved = clinicalPhase == 4) %>% 
  drop_na(Disease_class) %>% # Drops the few ontologies which couldn't be mapped
  mutate(in_OTG = symbol %in% otar.ibd.colocs$Molecular_trait) %>% 
  # Manually curated list of known drug targets
  filter(!symbol %in% c('ITGAL', 'ITGA4', 'JAK2', 'STAT3', 'TNFSF15','TNFRSF14'),
         Disease_class != 'Inflammatory bowel disease')

p <- ggplot(plot.summarised.drug.data, aes(x=mechanismOfAction, y=Disease_class, colour=in_OTG)) +
  # geom_image(aes(image=phase_image),size=0.04) +
  geom_point(shape=1,size=4,stroke=1.1) + 
  geom_point(data=filter(plot.summarised.drug.data, is_approved == TRUE),shape=19,size=2,
             show.legend = FALSE) +
  theme_classic() +
  scale_colour_manual(values=otar.palette,
                    labels=c('Novel','In OTG'),
                    name='Coloc novelty') +
  labs(x = "Mechanism of action", y = 'Disease group') +
  facet_grid(~symbol,scales='free',space='free') +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        legend.position = 'right',
        strip.clip = "off")
print(p)
if (SAVE.PLOTS) {
  ggsave(paste0(out.dir, 'doe_drugs.pdf'),
         plot = p, width = 21, height = 9,,device = cairo_pdf)
}
