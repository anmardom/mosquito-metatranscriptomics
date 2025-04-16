library(tidyr)
library(tidyverse)
setwd("/home/alumno/Descargas")
#fichero1 <- "/home/alumno/t/enrich/old_dataset5_joint_normalised_genefamilies.tsv"
fichero1 <- "old_dataset5_joint_normalised_genefamilies.tsv"
#fichero2 <- "/home/alumno/t/enrich/old_dataset5_joint_normalised_pathabundance.tsv"
fichero2 <- "old_dataset5_joint_normalised_pathabundance.tsv"

genefamilies <- read.delim(fichero1, sep = "\t")
#genefamilies$X..Gene.Family <- factor(genefamilies$X..Gene.Family)
pathabundance <- read.delim(fichero2, sep = "\t")
#pathabundance$X..Pathway <- factor(pathabundance$X..Pathway)

# Long format for pathabundances
datalong <- pathabundance %>%
  pivot_longer(cols = -1, names_to = "Muestra", values_to = "Abundancia") %>%
  mutate(grupo = case_when(
    grepl("SRR12064740|SRR12064741", Muestra) ~ "Midguts", 
    grepl("SRR12064742|SRR12064743", Muestra) ~ "Saliva"
  ))

# For 2 groups datasets
t_res <- datalong %>% 
  group_by(X..Pathway) %>% 
  summarise(
    p_value = tryCatch(
      t.test(Abundancia ~ grupo)$p.value, 
      error = function(e) NA
    )
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adj)

# For more than 2 groups datasets
anova_res <- datalong %>%
  group_by(X..Pathway) %>%
  summarise(
    p_anova = tryCatch(
      summary(aov(Abundancia ~ grupo))[[1]][["Pr(>F)"]][1], 
      error = function(e) NA
    )
  )

# Fold change for abundance comparisons
log2fc <- datalong %>%
  group_by(X..Pathway, grupo) %>%
  summarise(mean_abundancia = mean(Abundancia), .groups = "drop") %>%
  pivot_wider(names_from = grupo, values_from = mean_abundancia) %>%
  mutate(
    log2FC = log2(Midguts + 1) - log2(Saliva + 1)
  )

# Results joint
final_res <- left_join(t_res, log2fc, by = "X..Pathway") %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adj)

# Significative results
significativos <- final_res %>%
  filter(p_value < 0.05)

# Save table
write.table(final_res, file = "humann_results_ttest.tsv", sep = "\t", row.names = FALSE)

