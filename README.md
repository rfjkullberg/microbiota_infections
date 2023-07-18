# Gut microbiota and infectious disease hospitalisation

This is the code used for the main analyses in "Impact of butyrate-producing gut microbiota on the risk of infectious disease hospitalisation: results from two large population-based cohorts" (submitted). For questions: Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl

A file containing mock clinical metadata is included to show the type of data, data structure and its definitions ('Mock data and definitions.xlsx'). 
Data protection regulations do not allow public sharing of individual participant data on hospitalisations. Please refer to the data availability statement in the manuscript for information on how to gain data access. 

[![DOI](https://zenodo.org/badge/657511587.svg)](https://zenodo.org/badge/latestdoi/657511587)


## Step 1 - Load libraries

```
library(tidyverse)
library(data.table)
library(yingtools2)
library(phyloseq)
library(rlang)
library(microbiome)
library(vegan)
library(biomformat)
library(survival)
library(survminer)
library(cmprsk)
library(glmnet)
library(RColorBrewer) 
library(cowplot)
library(ggpubr) 
library(ggrepel)
library(scales)
library(readxl)
library(DESeq2)
library(ANCOMBC)
library(data.table)
library(MatchIt)
library(nlme)
library(rms)
```

## Step 2 - Load data

Microbiota sequence data is already preprocessed and a count table is produced, which is integrated with the taxonomy and a phylogenetic tree using the phyloseq package (details described in the manuscript). 

```
df <- read_csv("~/Documents/Helius/Data/metadata.csv") # metadata
P <- readRDS("~/Documents/PhD/Helius/Data/phyloseq.RDS") # phyloseq file 
```

We have defined 3 key variables: infection, event, end_date and time (see 'Data definitions' in the mock data file):
* Infection is a binary variable: participants with an admission for infection get "1"; those with mortality due to infection also get "1"; all others get "0". 
* Event is a variable with three levels: participants with admission/mortality due to infection get "infection"; those with mortality for other causes (competing risk) get "non-infection mortality"; all others get "no event"
* Time is the number of days between sample collection (start of follow-up) and the end date: The end date is the end of follow up: 31 December 2020 for participants without event (in the derivation cohort), the admission date (admission_date) for participants with hospital admission, and the mortality date (mortality_date) for participants with non-infection mortality (competing risk). Of note: mortality_date is the end of follow-up for participants with mortality due to infection, but without hospital admission.  

```
table(df$infection)
```

```
##    0    1 
## 4096  152
```

```
table(df$event)
```

```
##                no event               infection non-infection mortality 
##                    4046                     152                      50 
```

In steps 3-7, we assessed associations between our outcome variable and predefined key features of the bacterial gut microbiome: α-diversity, the relative abundance of butyrate-producing bacteria, and community composition. 


## Step 3 - Shannon diversity
Calculate alpha diversity

```
P.alpha <- microbiome::aggregate_taxa(P, level = "Species")
alpha <- estimate_richness(P.alpha)
alpha$sample_id <- row.names(alpha)                     
alpha$sample_id <- gsub("X","",as.character(alpha$sample_id))
alpha$sample_id <- str_replace(alpha$sample_id, "\\.","-")     
alpha <- alpha %>% dplyr::select(sample_id, Shannon)
df <- left_join(df, alpha) %>%
  mutate(tertiles = ntile(`Shannon`, 3)) %>%
  mutate(tertiles = as.factor(if_else(tertiles == 1, 'Low diversity', if_else(tertiles == 2, 'Intermediate diversity', 'High diversity')))) %>%
  mutate(tertiles = fct_relevel(tertiles, "Low diversity", "Intermediate diversity", "High diversity"))                
```

Dot plot showing distribution of Shannon diversity:
```
df %>%
  ggplot(aes(x= "", y = Shannon, fill = tertiles)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.2)+
  geom_hline(linetype = 2, max(df$Shannon[df$tertiles == "Intermediate diversity"])) + 
  geom_hline(linetype = 2, max(df$Shannon[df$tertiles == "Low diversity"])) + 
  scale_fill_manual(values=c( "#882255", "#DDCC77", '#44AA99')) +
  theme_cowplot()+
  xlab("")+
  ylab("Shannon diversity")+
  theme(legend.position = "none")                              
```

Competing risk analysis (cause-specific hazard ratio):
```
summary(coxph(Surv(time, event == "infection") ~ Shannon, data=df, id=sample_id)) # as continuous variable
summary(coxph(Surv(time, event == "infection") ~ tertiles, data=df, id=sample_id)) # in tertiles

table(df$tertiles, df$event)
```

Multivariable competing risk analysis (cause-specific hazard ratio):
```
summary(coxph(Surv(time, event == "infection") ~ 
                Shannon + age + sex + ethnicity +
                smoking + alcohol + prior_antibiotics +  physical_activity +
                comorb_diabetes + comorb_cardiovascular + comorb_cancer	+ comorb_hypertension + comorb_pulmonary + comorb_gastrointestinal, 
              data=df, id=sample_id))

summary(coxph(Surv(time, event == "infection") ~ 
                tertiles + age + sex + ethnicity +
                smoking + alcohol + prior_antibiotics +  physical_activity +
                comorb_diabetes + comorb_cardiovascular + comorb_cancer	+ comorb_hypertension + comorb_pulmonary + comorb_gastrointestinal, 
              data=df, id=sample_id))
```

We used the survminer package to plot the cumulative events per tertile of Shannon diversity
```
fit <- survfit(Surv(time, event == "infection") ~ tertiles, data=df)
ggsurvplot(fit, data = df, 
           fun = "event", 
           break.time.by = 365.25,  
           xscale = 365.25,
           axes.offset = F, 
           censor = F,            
           xlab = "Years since sample collection", 
           ylab = "Cumulative Incidence", 
           palette=c( "#882255", "#DDCC77", '#44AA99'))
```
```
rm(P.alpha, alpha, fit)
```

## 4 - Butyrate-producers 
We quantified the relative abundance of butyrate-producing bacteria by calculating the relative abundance of 15 bacteria that are known to be the most abundant drivers of butyrate production. This  might take a while to run.

```
butyrate.producers <- c("Butyricimonas", "Odoribacter", "Alistipes", "Eubacterium", "Anaerostipes","Butyrivibrio","Coprococcus", "Roseburia","Shuttleworthia","Butyricicoccus","Faecalibacterium","Flavonifractor","Pseudoflavonifractor","Oscillibacter", "Subdoligranulum")
butyrate.producers.species <- c("Subdoligranulum variabile")


butyrate<- get.otu.melt(P, filter.zero = F)
butyrate$Genus <- gsub("g__","",as.character(butyrate$Genus)) # remove the "g__" addition if necessary
butyrate$Species <- gsub("s__","",as.character(butyrate$Species)) # remove the "s__" addition if necessary


butyrate <- butyrate %>%
  splitstackshape::concat.split(split.col = "Genus", sep = "_") %>% # probably need to remove the additional Genus information (suffixes)
  mutate(Genus = Genus_1) %>%
  subset(Genus %in% butyrate.producers | Species %in% butyrate.producers.species) %>%
  group_by(sample_id)%>%                       
  summarize(sum = sum(pctseqs))
  
df <- left_join(df, butyrate) %>%
  mutate(butyrate = sum*100) %>% # convert fractions to percents
  mutate(tertiles_butyrate = ntile(`sum`, 3)) %>%
  mutate(tertiles_butyrate = as.factor(if_else(tertiles_butyrate == 1, 'Low butyrate', if_else(tertiles_butyrate == 2, 'Intermediate butyrate', 'High butyrate')))) %>%
  mutate(tertiles_butyrate = fct_relevel(tertiles_butyrate, "Low butyrate", "Intermediate butyrate", "High butyrate")) %>%
  dplyr::select(-sum)
```

Next, we generated a dot plot that shows the distribution of butyrate-producing bacteria:
```
df %>%
  ggplot(aes(x= "", y = butyrate, fill = tertiles_butyrate)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.2)+
  geom_hline(linetype = 2, max(df$butyrate[df$tertiles_butyrate == "Intermediate butyrate"])) + 
  geom_hline(linetype = 2, max(df$butyrate[df$tertiles_butyrate == "Low butyrate"])) + 
  scale_fill_manual(values=c( "#882255", "#DDCC77", '#44AA99')) +
  theme_cowplot()+
  xlab("")+
  ylab("Butyrate-producing bacteria")+
  theme(legend.position = "none")
```

Competing risk analysis (cause-specific hazard ratio):
```
summary(coxph(Surv(time, event == "infection") ~ butyrate, data=df, id=sample_id)) # as continuous variable
summary(coxph(Surv(time, event == "infection") ~ tertiles_butyrate, data=df, id=sample_id)) # in tertiles

table(df$tertiles_butyrate, df$event)
```

Multivariable competing risk analysis (cause-specific hazard ratio):
```
summary(coxph(Surv(time, event == "infection") ~ 
                butyrate + age + sex + ethnicity +
                smoking + alcohol + prior_antibiotics +  physical_activity +
                comorb_diabetes + comorb_cardiovascular + comorb_cancer	+ comorb_hypertension + comorb_pulmonary + comorb_gastrointestinal, 
              data=df, id=sample_id))

summary(coxph(Surv(time, event == "infection") ~ 
                tertiles_butyrate + age + sex + ethnicity +
                smoking + alcohol + prior_antibiotics +  physical_activity +
                comorb_diabetes + comorb_cardiovascular + comorb_cancer	+ comorb_hypertension + comorb_pulmonary + comorb_gastrointestinal, 
              data=df, id=sample_id))
```

We used the survminer package to plot the cumulative events per tertile of butyrate-producers:
```
fit <- survfit(Surv(time, event == "infection") ~ tertiles_butyrate, data=df)
ggsurvplot(fit, data = df, 
           fun = "event", 
           break.time.by = 365.25, 
           xscale = 365.25,
           axes.offset = F, 
           censor = F,
           #xlim = c(0,2600),
           show.legend = F,
           xlab = "Years since sample collection", 
           ylab = "Cumulative Incidence", 
           palette = c( "#882255", "#DDCC77",  "#44AA99"))
```

We performed a CLR-transformation to correct for the compositional nature of microbiome data
```
P.clr <- P
tax.CLR <- get.tax(P.CLR) 
tax.CLR$Genus <- gsub("g__","",as.character(tax.CLR$Genus)) # remove the "g__" addition if necessary
tax.CLR$Species <- gsub("s__","",as.character(tax.CLR$Species)) # remove the "s__" addition if necessary

# Sum up the counts for each of the butyrate-producers and the rest to "Other
tax.CLR <- tax.CLR %>%
  splitstackshape::concat.split(split.col = "Genus", sep = "_") %>% # probably need to remove the additional Genus information
  mutate(Genus = as.character(Genus_1))%>%
  mutate(Species = as.character(Species)) %>%    
  mutate(butyrateproducers = if_else(Genus %in% butyrate.producers, Genus, 
                               if_else(Species %in% butyrate.producers.species, Species, "Other"))) %>%
  dplyr::select(otu, butyrateproducers)
tax_table(P.CLR) <- set.tax(tax.CLR)
P.CLR <- microbiome::aggregate_taxa(P.CLR, level = "butyrateproducers")

# CLR-transformation
P.CLR <- microbiome::transform(P.CLR, "clr")

butyrate.clr <- get.otu.melt(P.CLR, filter.zero = F) %>%
  filter(butyrateproducers != "Other") %>%
  group_by(sample)%>%
  summarize(butyrate.clr = sum(numseqs)) %>% # calculate the CLR-transformed abundance of the butyrate-producing bacteria per sample
  mutate(tertiles_butyrate_clr = ntile(`butyrate.clr`, 3)) %>%
  mutate(tertiles_butyrate_clr = as.factor(if_else(tertiles_butyrate_clr == 1, 'Low butyrate', if_else(tertiles_butyrate_clr == 2, 'Intermediate butyrate', 'High butyrate')))) %>%
  mutate(tertiles_butyrate_clr = fct_relevel(tertiles_butyrate_clr, "Low butyrate", "Intermediate butyrate", "High butyrate")) %>%
  dplyr::rename(sample_id=sample)

df <- left_join(df, butyrate.clr)
```

```
table(df$event, df$tertiles_butyrate_clr)
summary(coxph(Surv(time, event == "infection") ~ butyrate.clr, data=df, id=sample_id)) # as continuous variable
summary(coxph(Surv(time, event == "infection") ~ tertiles_butyrate_clr, data=df, id=sample_id)) # in tertiles
```
```
summary(coxph(Surv(time, event == "infection") ~ 
                butyrate.clr + age + sex + ethnicity +
                smoking + alcohol + prior_antibiotics +  physical_activity +
                comorb_diabetes + comorb_cardiovascular + comorb_cancer	+ comorb_hypertension + comorb_pulmonary + comorb_gastrointestinal, 
              data=df, id=sample_id))

summary(coxph(Surv(time, event == "infection") ~ 
                tertiles_butyrate_clr + age + sex + ethnicity +
                smoking + alcohol + prior_antibiotics +  physical_activity +
                comorb_diabetes + comorb_cardiovascular + comorb_cancer	+ comorb_hypertension + comorb_pulmonary + comorb_gastrointestinal, 
              data=df, id=sample_id))
```

```
rm(butyrate, butyrate.clr, otu, OTU, P.CLR, tax.CLR, butyrate.producers, butyrate.producers.species, fit)
```

## 5 - Beta diversity
```
P.comp <- microbiome::transform(P, "compositional")
set.seed(88)
bray <- phyloseq::distance(P.comp, method = "bray") 
adonis2((bray ~ infection), by = 'margin', data = df, permutations =999)
```

We created a figure showing the centroids of our two main groups (infection / no infection) and a figure coloured by the abundance of butyrate-producing bacteria.
```
ord <- ordinate(P.comp, method = "PCoA", distance = bray)

# Create dataframe with coordinates and centroids of participants
bray_df <- plot_ordination(P.comp, ord, type = "samples", color = "infection", justDF = T) %>%
  dplyr::select(Axis.1, Axis.2, sample_id, infection)
centroids <- aggregate(cbind(Axis.1, Axis.2)~infection, data=bray_df, mean)
bray_df <- merge(bray_df, centroids, by="infection", suffixes=c("", ".centroid"))

# Create separate dataframes for participants with/without infection
bray_df_inf <- bray_df %>%
  filter(infection == "1")   
bray_df_no <- bray_df %>%   
  filter(infection == "0")
  
ggplot() +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) +
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) +
  geom_point(data = bray_df, aes(x=Axis.1, y=Axis.2), size = 1, colour = "grey36") +
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, colour=infection, label=c("No infection", "Infection")), size=6, fill=ggplot2::alpha(c("white"),.76)) +
  stat_ellipse(data = bray_df_inf, aes(Axis.1, Axis.2, colour = infection), level = 0.95, lwd= 3) +
  stat_ellipse(data = bray_df_no, aes(Axis.1, Axis.2, colour = infection), level = 0.95, lwd= 3) +
  geom_point(data = centroids, aes(x=Axis.1, y=Axis.2, colour = infection), size = 5) +
  theme_cowplot() +
  scale_colour_manual(values=c("#4197CC", "#9F514D")) +
  xlab("Axis1")+
  ylab("Axis2")+
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank())
```
```
# Colour by butyrate
bray_df <- bray_df %>%
  left_join(df)

ggplot() +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) +
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) +
  geom_point(data = bray_df, aes(x=Axis.1, y=Axis.2, colour = butyrate), size = 1.5) +
  stat_ellipse(data = bray_df_inf, aes(Axis.1, Axis.2), level = 0.95, lwd= 3) +
  stat_ellipse(data = bray_df_no, aes(Axis.1, Axis.2), level = 0.95, lwd= 3) +
  geom_point(data = centroids, aes(x=Axis.1, y=Axis.2), size = 5) +
  theme_cowplot() +
  xlab("Axis1")+
  ylab("Axis2")+
  scale_colour_gradient2(low = "#821d4e", mid = "#DDCC77", high = "#1d8251", midpoint = median(bray_df$butyrate))+
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank())
```
```
rm(P.comp, bray, ord, deseq, centroids, bray_df, bray_df_inf, bray_df_no, bray_dsq)
```

## 6 - DESeq2
We used DESeq2 to identify differentially abundant genera between participants with an infection (or infection related mortality) during follow-up and those without infection.

Derivation (HELIUS) cohort only:
```
ps.dsq <- tax_glom(P, "Genus")
ps.dsq <- core(ps.dsq, detection=1, prevalence=10/100, include.lowest=T) # keep genera prevalent in at least 10% of the population

gm_mean <- function(x, na.rm=T){
  exp(sum(log(x[x>0]), na.rm=na.rm)/length(x))
}

dsq <- phyloseq_to_deseq2(ps.dsq, ~infection)
geoMeans <- apply(counts(dsq), 1, gm_mean)
dsq <- estimateSizeFactors(dsq, geoMeans = geoMeans) 
dsq <- DESeq(dsq,  fitType="local")  
res <- results(dsq, cooksCutoff = FALSE, pAdjustMethod = "BH" ) # Benjamini-Hochberg (BH) adjustment for multiple testing
deseq <- res[which(res$padj<0.05),] # select genera with BH-corrected P<0.05
deseq <- cbind(as(deseq, "data.frame"), as(tax_table(ps.dsq)[rownames(deseq),], "matrix"))

deseq <- deseq %>%
  group_by(Genus) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T) %>%
  mutate(group=ifelse(log2FoldChange < 0, "Infection", "No infection"))

ggplot(deseq, aes(x=reorder(Genus, log2FoldChange), y=log2FoldChange, fill=group)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  ylab("log 2-Fold change") +
  xlab("") + 
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("#4197CC", "#9F514D")) +
  theme(legend.position = "none")
```
Taxa that were differentially abundant (BH-adjusted p<0.05) in the HELIUS cohort, were validated in the FINRISK cohort. 
Validation (FINRISK) cohort only:
```
ps.dsq <- tax_glom(P, "Genus")
ps.dsq <- core(ps.dsq, detection=1, prevalence=10/100, include.lowest=T) # keep genera prevalent in at least 10% of the population

gm_mean <- function(x, na.rm=T){
  exp(sum(log(x[x>0]), na.rm=na.rm)/length(x))
}

dsq <- phyloseq_to_deseq2(ps.dsq, ~infection)
geoMeans <- apply(counts(dsq), 1, gm_mean)
dsq <- estimateSizeFactors(dsq, geoMeans = geoMeans) 
dsq <- DESeq(dsq,  fitType="local")  
res <- results(dsq, cooksCutoff = FALSE, pAdjustMethod = "BH" ) # Benjamini-Hochberg (BH) adjustment for multiple testing

# Genera with significant DESeq2 results in HELIUS:
helius.sign.deseq.gen <- c("Butyrivibrio_A_180067", "Streptococcus", "Acutalibacteraceae bacterium", "Veillonella_A",
                           "Merdicola", "G11", "CAG-307", "Stercorousia",
                           "COE1",   "CAG-1427", "Zag111", "Hydrogeniiclostridium", 
                           "Megamonas",  "Anaerobutyricum", "Ruminococcaceae bacterium")

deseq <- cbind(as(res, "data.frame"), as(tax_table(P.dsq)[rownames(res),], "matrix")) # results from DESeq2
deseq$Family <- gsub("f__", "", as.character(deseq$Family)) # remove the "f__" addition if necessary
deseq$Genus <- gsub("g__", "", as.character(deseq$Genus)) # remove the "g__" addition if necessary
deseq <- deseq %>%
  mutate(Genus = if_else(Genus == "", paste(Family, "bacterium", sep = " ", collapse = NULL), Genus)) %>%
  dplyr::filter(Genus %in% helius.sign.deseq.gen)
print(deseq)

deseq <- deseq %>%
  dplyr::filter(pvalue < 0.05) %>%
  group_by(Genus) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T) %>%
  mutate(group=ifelse(log2FoldChange < 0, "Infection", "No infection")) 

deseqplot <- ggplot(deseq, aes(x=reorder(Genus, log2FoldChange), y=log2FoldChange, fill=group)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  ylab("log 2-Fold change") +
  xlab("") + 
  theme_bw() +
  scale_fill_manual(values = c("#4197CC", "#9F514D")) +
  theme(legend.position = "none")
```

Statistical methods that aim to identify differentially abundant taxa in microbiota research are vulnerable to false positive findings [https://doi.org/10.1038/s41467-022-28034-z]. We therefore used Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC2) instead of DESeq2 to confirm differentially abundant genera. 

```
out <- ancombc2(ps.dsq, 
                fix_formula="infection", 
                tax_level = "Genus", 
                p_adj_method = "BH", 
                prv_cut = 0)

ancom <- out$res %>%
  mutate(lfc = lfc_infection1) %>%
  filter(q_infection1 < 0.05) %>%
  dplyr::select(taxon, lfc, q_infection1) %>%
  mutate(group = if_else(lfc < 0, "Infection", "No infection")) 

ggplot(ancom, aes(x=reorder(taxon, lfc), y=lfc, fill=group)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  ylab("Log fold change") +
  xlab("") + 
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("#9F514D", "#4197CC")) +
  theme(legend.position = "none")
```

## 7 - Risk score
Since DESeq2 and ANCOM-BC show differences in gut microbiota composition on the group level, rather than the effect of changes in an individual participant on his/her risk of infection, we sought to identify a signature of bacterial abundances associated with the risk of infection.
We used regularized Cox regression with cross-validation to derive a risk score as earlier described: Peled et al. NEJM 2020 [https://doi.org/10.1056/NEJMoa1900623]

The signature of bacterial effect sizes (coefficient of each term in the model) was trained in the HELIUS cohort using regularized regression. Positive values indicate increased risk of infection; negative values decreased risk. The risk score was subsequently calculated by multiplying the weights (defined by the regularized Cox model) with the relative abundances. The risk score trained in the HELIUS cohort was tested in the validation cohort of participants from the FINRISK study. 
A full list of bacteria contributing to the risk score (including their weights) are provided in 'risk-score-helius-infections.xlsx'. 

To calculate the risk score (HELIUS cohort only):
```
P.cox <- P
tax <- tax %>%
  mutate(Genus = if_else(Genus == "", paste(Family, "bacterium", sep = " ", collapse = NULL), Genus)) %>%
  filter(Genus != " bacterium")
tax_table(P.cox) <- set.tax(tax)
P.cox <- phyloseq::tax_glom(P.cox, "Genus", NArm=T)
P.cox <- core(P.cox, detection=1, prevalence=10/100, include.lowest=T)
P.cox <- microbiome::transform(P.cox, "compositional")
tax <- get.tax(P.cox)
  
# Define response and predictors variables
response <- Surv(df$time, df$event=="infection") 
predictors <- otu_table(P.cox)
predictors <- log(predictors+2e-05) # add a pseudo count to eliminate the possibility of –Infinity values

# Calculate risk score
set.seed(88)
fit <- glmnet::cv.glmnet(predictors, response, 
                         alpha=c(0, .1, .25, .5, 1.0), 
                         family="cox", 
                         maxit=10000)
                         
riskscore <- as.data.frame(as.matrix(coef(fit, s=fit$lambda.1se))) %>%
  rownames_to_column(var='otu') %>%
  mutate(weight = `1`) %>%
  left_join(tax) %>%
  dplyr::select(Phylum, Class, Order, Family, Genus, weight) %>%
  filter(Genus != " bacterium") %>%
  mutate_at(c("weight"), ~(scale(.))) %>%
  mutate(direction = if_else(weight > 0, "infection", "no infection"))              
```

Visualizing the risk score (HELIUS cohort only):
```
scaled_riskscore <- get.otu.melt(P.cox, filter.zero = F) %>%
  group_by(Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarize(meanabundance = mean(pctseqs)) %>%
  left_join(riskscore) %>%
  mutate(abundance_scaled_weight = meanabundance * abs(weight)) %>%
  mutate(Genus = as.character(Genus)) %>%
  mutate(label = if_else(abundance_scaled_weight>0.01, Genus, 
                         if_else(abs(weight)>2, Genus, ""))) %>%
  mutate(size = if_else(meanabundance > 0.05, "0.05", 
                        if_else(meanabundance > 0.01, ">0.01", 
                                if_else(meanabundance > 0.001, ">0.001", "<0.001"))))

ggplot(scaled_riskscore, aes(x=weight, y=abundance_scaled_weight, colour=direction)) + 
  geom_point(aes(size=size)) +
  geom_text_repel(aes(label=label), min.segment.length=unit(0, "lines"), 
                  nudge_x=0, nudge_y=0.2, segment.alpha=0.3, force=2) +
  scale_y_continuous(trans=log_epsilon_trans(epsilon=.01)) +
  scale_colour_manual(values = c("#9F514D", "#4197CC")) +
  theme_cowplot() +
  theme(legend.position = "none")
```

To validate the risk score (FINRISK cohort only), we first calculated the risk score for each individual: 
```
riskscores <- read_excel("risk-score-helius-infections.xlsx") # file to calculate the risk scores

riskscore <- get.otu.melt(P)
riskscore <- riskscore %>%
  mutate(Genus = if_else(Genus == "", paste(Family, "bacterium", sep = " ", collapse = NULL), Genus)) %>%
  filter(Genus != " bacterium") %>%
  left_join(riskscores, by = Genus) %>%  
  filter(!is.na(weight)) %>% 
  mutate(indiv_score = weight*pctseqs) %>%
  group_by(sample_id, infection) %>%                                      
  summarize(indiv_score = sum(indiv_score)) %>% 
  left_join(df) %>% 
  ungroup() %>%
  mutate(tertiles_riskscore = ntile(indiv_score, 3)) %>%
  mutate(tertiles_riskscore = as.factor(if_else(tertiles_riskscore == 1, 'Low risk score', if_else(tertiles_riskscore == 2, 'Intermediate risk score', 'High risk score')))) %>%
  mutate(tertiles_riskscore = fct_relevel(tertiles_riskscore, "Low risk score", "Intermediate risk score", "High risk score")) 
```

Dot plot showing distribution of risk scores:
```
riskscore %>%
  ggplot(aes(x= "", y = indiv_score, fill = tertiles_riskscore)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.2)+
  geom_hline(linetype = 2, max(riskscore$indiv_score[riskscore$tertiles_riskscore == "Intermediate risk score"])) + 
  geom_hline(linetype = 2, max(riskscore$indiv_score[riskscore$tertiles_riskscore == "Low risk score"])) + 
  scale_fill_manual(values=c("#44AA99", "#DDCC77","#882255")) +
  theme_cowplot()+
  xlab("")+
  ylab("Risk score")+
  theme(legend.position = "none")
```

Competing risk analysis (cause-specific hazard ratio):
```
summary(coxph(Surv(time, event == "infection") ~ indiv_score, data=df, id=sample_id)) # as continuous variable
summary(coxph(Surv(time, event == "infection") ~ tertiles_riskscore, data=df, id=sample_id)) # in tertiles

table(df$tertiles, df$event)
```

We used the survminer package to plot the cumulative events per tertile of risk score
```
fit <- survfit(Surv(time, event == "infection") ~ tertiles_riskscore, data=riskscore)
survp <- ggsurvplot(fit, data = riskscore, 
           fun = "event", 
           break.time.by = 365.25,    
           xscale = 365.25,
           axes.offset = F, 
           censor = F, 
           show.legend = F,
           #xlim = c(0,2600),
           xlab = "Years since sample collection", 
           ylab = "Cumulative Incidence", 
           palette = c("#44AA99", "#DDCC77","#882255"))
```

## 8 - Matched analyses
Nested age-, sex-, ethnicity-, antibiotics-, and comorbidity-matched case-control analyses within the derivation cohort were used to examine differences in gut microbiota between participants with an infection-related hospitalisation (cases) and matched controls.

```
m <- df %>%
  mutate(age_groups = cut(age, c(18, seq(40, 60, by=10), Inf), include.lowest=T)) %>%
  filter(prior_antibiotics != "NA") %>%
  filter(comorb_pulmonary != "NA") %>%
  filter(comorb_gastrointestinal != "NA") %>%
  filter(comorb_diabetes != "NA")
  # 146 participants left; six participants with an infection had missing data on matching factors
```
```
set.seed(88)
matched <- matchit(infection ~ sex+ethnicity+age_groups+prior_antibiotics+comorb_pulmonary+comorb_gastrointestinal+comorb_diabetes, 
                    data = m, 
                    method = "nearest", 
                    ratio=1)
                    
matched <- match.data(matched)
```

First, we compared Shannon diversity and the abundance of butyrate-producing bacteria between cases and matched controls. 
```
P.matched <- P
sample_data(P.matched) <- set.samp(matched)

alpha <- estimate_richness(P.matched)
alpha$sample <- row.names(alpha)
alpha$sample <- gsub("X","",as.character(alpha$sample))
alpha <- alpha %>% select(sample, Shannon)
matched <- left_join(matched, alpha) 

comparisons <- list(c("1", "0")) 
matched %>%
  ggplot(aes(x = infection, y = Shannon, fill = infection)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.0) +
  theme_bw() +
  xlab("") +
  ylab("Shannon diversity") +
  scale_fill_manual(values = c("#4197cc", "#9f514d")) +
  theme(legend.position = "none")
anova(lme(Shannon ~ infection, random=~1|subclass,  data = matched))
```

```
matched %>%
  ggplot(aes(x = infection, y = butyrate, fill = infection)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.0) +
  theme_bw() +
  xlab("") +
  ylab("Butyrate-producers") +
  scale_fill_manual(values = c("#4197cc", "#9f514d")) +  
  theme(legend.position = "none")
anova(lme(butyrate ~ infection, random=~1|subclass,  data = matched))
```

Next, we tested for differences in community composition

```
P.comp <- microbiome::transform(P.matched, "compositional")
set.seed(88)
bray <- phyloseq::distance(P.comp, method = "bray")
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(matched)
adonis2((bray ~ infection), by = 'margin', data = df.bray, permutations =999)
ord <- ordinate(P.comp, method = "PCoA", distance = bray) 

# Create dataframe with coordinates and centroids of participants 
bray_df <- plot_ordination(P.comp, ord, type = "samples", color = "infection", justDF = T) 
centroids <- aggregate(cbind(Axis.1, Axis.2)~infection, data=bray_df, mean)
bray_df <- merge(bray_df, centroids, by="infection", suffixes=c("", ".centroid"))

# Create separate dataframes for participants with/without infection 
bray_df_inf <- bray_df %>%
  filter(infection == "1")
bray_df_no <- bray_df %>%
  filter(infection == "0")

ggplot() +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) +
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) +
  stat_ellipse(data = bray_df_inf, aes(Axis.1, Axis.2, colour = infection), level = 0.95, lwd= 3) +
  stat_ellipse(data = bray_df_no, aes(Axis.1, Axis.2, colour = infection), level = 0.95, lwd= 3) +
  geom_point(data = centroids, aes(x=Axis.1, y=Axis.2, colour = infection), size = 5) + geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, colour=infection, label=c("No infection", "Infection")), size=6, fill=ggplot2::alpha(c("white"),.76)) +
  theme_cowplot() +
  scale_colour_manual(values=c("#4197CC", "#9F514D")) +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.text = element_blank())
```

## 9 - Role of covariates
In exploratory analyses in the derivation cohort, we further studied potential associations between covariates (such as demographics, lifestyle factors, antibiotics and comorbidities), the abundance of butyrate-producing bacteria, and infections. 

First, we explored associations between covariates (sex and age are shown here as an example) and the abundance of butyrate-producing bacteria. 
```
comparisons <- list(c("1", "2"))
df %>%
  ggplot(aes(x= sex, y = butyrate, fill = sex)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 2, width=0.15) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.format", size=5)+
  scale_fill_manual(values=c('#88CCEE', "#cc6677")) +
  theme_bw()+
  xlab("")+
  ylab("Butyrate-producing bacteria (%)")+
  ggtitle("Sex")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
  
df %>%
  ggplot(aes(x = age, y = butyrate))+
  geom_point(size=2,  colour = "#8B8989") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw()+
  xlab("Age")+
  ylab("Butyrate-producing bacteria (%)")+
  ggtitle("Age")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
summary(lm(age~butyrate, data=df)) 
```
Similar code was used to assess the other covariates (ethnicity, BMI, alcohol usage, smoking, physical activity, antibiotic exposure, comorbidities, and dietary variables). 

Next, we assessed univariable associations between these covariates and the risk of infection, using competing risk regression models. 
```
surv_object <- Surv(df$time, df$event == "infection")

cshr <- as.data.frame(rbind(
      summary(coxph(surv_object ~ age, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ sex, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ ethnicity, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ bmi_groups, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ smoking, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ alcohol, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ physical_activity, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ prior_antibiotics, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ comorb_diabetes, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ comorb_cardiovascular, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ comorb_cancer, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ comorb_hypertension, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ comorb_pulmonary, data = df, id = sample_id))[["conf.int"]], 
      summary(coxph(surv_object ~ comorb_gastrointestinal, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ diet_fattyacids, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ diet_satfattyacids, data = df, id = sample_id))[["conf.int"]],
      summary(coxph(surv_object ~ diet_fibres, data = df, id = sample_id))[["conf.int"]]))

pval <- as.data.frame(rbind(
  summary(coxph(surv_object ~ age, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ sex, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ ethnicity, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ bmi_groups, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ smoking, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ alcohol, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ physical_activity, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ prior_antibiotics, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ comorb_diabetes, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ comorb_cardiovascular, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ comorb_cancer, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ comorb_hypertension, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ comorb_pulmonary, data = df, id = sample_id))[["coefficients"]], 
  summary(coxph(surv_object ~ comorb_gastrointestinal, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ diet_fattyacids, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ diet_satfattyacids, data = df, id = sample_id))[["coefficients"]],
  summary(coxph(surv_object ~ diet_fibres, data = df, id = sample_id))[["coefficients"]])) %>%
  mutate(pvalue = as.numeric(`Pr(>|z|)`)) %>%
  dplyr::select(pvalue)
  
# combine
coxmodels <- as.data.frame(cbind(cshr, pval)) %>%
  mutate(cshr = as.numeric(`exp(coef)`)) %>%
  mutate(lower = as.numeric(`lower .95`)) %>%
  mutate(upper = as.numeric(`upper .95`)) %>%
  rownames_to_column(var = "variable") %>%
  dplyr::select(variable, cshr, lower, upper, pvalue) %>%
  rbind(list(variable = "Dutch ethnicity", cshr=1.00, lower=1.00, upper=1.00, pvalue=1.00)) %>% # setting 'Dutch ethnicity' as reference value
  rbind(list(variable = "Smoking, never", cshr=1.00, lower=1.00, upper=1.00, pvalue=1.00)) %>% # setting 'Smoking, never' as reference value
  mutate(variable = fct_relevel(variable, 
                                "diet_fibres", "diet_satfattyacids", "diet_fattyacids",
                                "comorb_cancer", "comorb_gastrointestinal1", "comorb_pulmonary1", "comorb_hypertension1", "comorb_cardiovascular1", "comorb_diabetes1", 
                                "med_corticosteroids1", "probiotics1", "prior_antibiotics1", "physical_activity1", "alcohol2", "smoking3","smoking1", "Smoking, never",
                                "bmi_groups≥30", "bmi_groups25-29", "BMI <25"
                                "ethnicity5", "ethnicity4", "ethnicity3","ethnicity2_2", "ethnicity2_1", "Dutch ethnicity",
                                "sex2", "age" ))

coxmodels %>%
  ggplot(aes(x = variable, y = cshr, ymin = lower, ymax = upper)) +
  geom_pointrange() + 
  geom_text(aes(y = -.2, x = variable, label = paste("p=", round(pvalue, 4)))) +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  xlab("") +
  ylab("Cause-specific hazard ratio (95% confidence interval)") +
  theme_bw()
```

```
rm(surv_object, cshr, pval coxmodels)
```

## 10 - Contrasts
We computed contrasts to assess differences between subgroups in the association between infection risk and abundances of butyrate-producing bacteria (i.e. whether covariates modified the relation between butyrate-producers and our primary outcome).

```
dd <- datadist(df)
options(datadist='dd')

summary(coxph(Surv(time, event=="infection") ~ Shannon * sex, df))
cox.Shannonsex <- cph(Surv(time, event == "infection") ~ Shannon * sex, data=df)
p.Shannonsex <- as.data.frame(Predict(cox.Shannonsex, Shannon, sex, fun = exp, ref.zero = T))

ggplot(p.Shannonsex, aes(x=Shannon, colour = sex))+
  geom_line(aes(y = yhat), size=2) +
  scale_y_continuous(name = "Risk of infection (cause-specific Hazard Ratio)",
                     trans = "log")) +
  scale_x_continuous(name = "") +
  scale_colour_manual(values=c('#88CCEE', "#cc6677")) +
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=(1), linetype="dashed", color = "grey36", size=1)

cit_crude<-cph(formula=Surv(time,event == "infection") ~ 
                 Shannon*sex, data=df,x=T, y=T)
kit <- seq(2.5,5,by=0.25) 
w <- contrast(cit_crude, 
              list(sex='2', Shannon=kit), 
              list(sex='1', Shannon=kit))
w <- cbind(kit, w[["Contrast"]], w[["Lower"]], w[["Upper"]], w[["Pvalue"]])
colnames(w) <- c("Shannon", "Contrast", "Lower", "Upper", "Pvalue")
w <- as.data.frame(w)

w %>%
  ggplot(aes(x=Shannon, y=exp(Contrast))) +
  geom_pointrange(aes(ymin=exp(Lower), ymax=exp(Upper)), width=.2) +
  annotate(geom="text",  x=Inf, y = Inf, hjust=1, vjust=1.1, size=4,colour = "black",
           label=c(paste("Risk for females relative to males ")))+
  geom_hline(yintercept = 1) +
  scale_y_continuous(name = "",
                     trans = "log") +
  theme_bw()
```
Similar code was used to assess the other covariates. 

