# Rscript metagenomeSeq differential abundance analysis
# Angel Martin Dominguez

# Load libraries
library(metagenomeSeq)
library(ggplot2)

# Load datasets files required for MRexperiment object
# Dataset 1: Plasmodium infection
matrix <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset1_Exp2_counts.tab"
OTUS <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset1_Exp2_OTU.tsv"
metadata <- "/home/alumno/Descargas/TFM_Angel/res2/metadata1_infection_control.tsv"

# Dataset 2: Serratia marcescens resistance
matrix <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset2_Exp4_counts.tab"
OTUS <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset2_Exp4_OTU.tsv"
metadata <- "/home/alumno/Descargas/TFM_Angel/res2/metadata2_single_serratia_control.tsv" # Metadata 1
metadata <- "/home/alumno/Descargas/TFM_Angel/res2/metadata2_serratias_control.tsv"       # Alternative metadata 2

# Dataset 3: Trypanosoma brucei susceptibility
matrix <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset3_Exp5_counts.tab"
OTUS <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset3_Exp5_OTU.tsv"
metadata <- "/home/alumno/Descargas/TFM_Angel/res2/metadata3_single_tryp_control.tsv"

# Dataset 4: Insecticide resistance
matrix <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset4_Exp3_counts.tab"
OTUS <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset4_Exp3_OTU.tsv"
metadata <- "/home/alumno/Descargas/TFM_Angel/res2/metadata4_bus-sus-res_kis.tsv"

# Dataset 5: Infected tissues comparison
matrix <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset5_Exp1_counts.tab"
OTUS <- "/home/alumno/Descargas/TFM_Angel/res2/Dataset5_Exp1_OTU.tsv"
metadata <- "/home/alumno/Descargas/TFM_Angel/res2/metadata5_midguts_salivary.tsv"


# Load and prepare data for phys object
sample_matrix <- read.delim(matrix, sep = "\t", row.names = 1)    # Read count matrix
otus <- read.delim(OTUS, sep = "\t", row.names = 1)               # Read otu file
rownames(sample_matrix) <- otus$OTU                               # Replace rownames with OTU names in count matrix
rownames(otus) <- otus$OTU                                        # Replace rownames with OTU names in otu table
mdt <- read.delim(metadata, sep = "\t", row.names = 1)            # Read metadata
#mdt <- mdt[-8,]                                                  # ONLY in Dataset 3 (Exp 5), Discard inf_24 outlier
mdt$Concentration <- colSums(sample_matrix)                       # Replace metadata names with sample filenames
rownames(mdt) <- colnames(sample_matrix)
#sample_matrix <- sample_matrix[, -c(15, 16)]                     # ONLY in Dataset 4, Delete duplicated columns count matrix
#mdt <- mdt[-c(15, 16), ]                                         # ONLY in Dataset 4, Delete duplicated columns metadata

# Prepare elements and build MRexperiment object for analysis
pdata <- AnnotatedDataFrame(mdt)
otu_tabla <- AnnotatedDataFrame(otus)
mrexp <- newMRexperiment(sample_matrix, phenoData = pdata, featureData = otu_tabla)

# Inspect library size of samples and graph them in increasing size order
x <- as.data.frame(pData(mrexp))
x <- x[order(x$Concentration),]
x$Index <- seq(nrow(x))
ggplot(x, aes(x = Index, y = Concentration, color = Source)) + geom_point()

# Filter rare species that are present in few replicas
nrow(MRcounts(mrexp))
Num <- min(table(pData(mrexp)$Source))                        # Num is the number of replicas from lowest condition count
rareFeatures <- which(rowSums(MRcounts(mrexp) > 0) <= Num)
length(rareFeatures)
rF <- featureData(mrexp[rareFeatures, ])$OTU
# Save rareFeatures to check discarded species
#writeLines(rF, "/home/alumno/Descargas/TFM_Angel/res2/final/rareFeatures/dataset2_full_rareFeatures.txt")

# Delete rare species from mrexp
if(length(rareFeatures) > 0){
  mrexp = mrexp[-rareFeatures, ]
}

# Normalization of mrexp with scaling factors and proper quantile (in most cases 0.5)
normFactors(mrexp) <- rnorm(ncol(mrexp))
p <- cumNormStat(mrexp)
mrexp_norm <- cumNorm(mrexp, p = p)

# PCoA normalized counts
ncounts <- MRcounts(mrexp_norm, norm = TRUE)
log_ncounts <- log2(ncounts + 1)
pca_res <- prcomp(t(log_ncounts), center = TRUE, scale. = TRUE)
varianza <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
Source <- as.factor(pData(mrexp_norm)$Source)
pca_data <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Condition = Source)

# Create PCoA object
pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) + 
  geom_point(size = 3) + 
  labs(x = paste0("PC1 (", round(varianza[1], 1), "% variance)"), 
       y = paste0("PC2 (", round(varianza[2], 1), "% variance)"), 
       title = "PCA log2 normalized counts") + 
  theme_minimal()
print(pca)

# Save PCoA
#png("/home/alumno/Descargas/TFM_Angel/res2/final/pca/pca_dataset2_rep4_full_serratiaVScontrol.png", width = 5000, height = 6000, res = 600)
#print(pca)
#dev.off()


# Factors for differential abundance model matrix
normFactor <- normFactors(mrexp_norm)
normFactor <- log2(normFactor/median(normFactor) + 1)
Source<-as.factor(pData(mrexp_norm)$Source)

# fitZIG (fit for Zero-inflated Gaussian model)
modelo <- model.matrix(~0+Source+normFactor)      # Model matrix for fitZIG with factors Source and normFactor
ajustes = zigControl(maxit = 10, verbose = F)     # Model settings, if fzm fails, change maxit value to 1
fzm <- fitZig(obj = mrexp_norm, mod = modelo, useCSSoffset = FALSE, control = ajustes)  # fitZIG

# Extract fit and design model information for contrast matrixes required for pairwise comparison of conditions
zigFit = slot(fzm, "fit")
finalMod = slot(fzm, "fit")$design

# Get number of tax in MRexp
tax <- sapply(strsplit(as.character(fData(mrexp)$OTU), split = ";"), 
              function(i) {
                i[length(i)]
              })
# Print coefficients of fitZIG model to get coefficient names for contrast matrixes
head(MRcoefs(fzm, taxa = tax, coef = c(1:3)))    # Length c1:4 in dataset 3


# Contrasts matrixes to extract comparisons between conditions

# Contrast Dataset 1: Plasmodium infection
contrast.matrix = makeContrasts(Sourceinfection - Sourcecontrol, levels = finalMod)

# Contrasts Dataset 2: Serratia marcescens resistance
# First metadata table: Serratia administration + Time of sampling
contrast.matrix = makeContrasts(a = (Sourceserratia_24 - Sourceserratia_0), b = (Sourceserratia_0 - Sourcecontrol_0),
                                c = (Sourceserratia_24 - Sourcecontrol_24), 
                                levels = finalMod)

# Second metadata table: Serratia administration
contrast.matrix = makeContrasts(Sourceserratia - Sourcecontrol, levels = finalMod)

# Contrasts Dataset 3: Trypanosoma brucei susceptibility
contrast.matrix = makeContrasts(a = (Sourceinfection_24 - Sourcecontrol_24), b = (Sourceinfection_48 - Sourcecontrol_48), 
                                c = (Sourceinfection_48 - Sourceinfection_24), d = (Sourcecontrol_48 - Sourcecontrol_24), 
                                levels = finalMod)

# Contrasts Dataset 4: Insecticide resistance
contrast.matrix = makeContrasts(a = (Sourcesusceptible - Sourceresistant), b = (Sourcesusceptible - Sourcecontrol), 
                                c = (Sourceresistant - Sourcecontrol), 
                                levels = finalMod)

# Contrast Dataset 5: Infected tissues comparison
contrast.matrix = makeContrasts(Sourcesalivary - Sourcemidguts, levels = finalMod)


# Contrast fit for fitZIG + eBayes correction
fit = contrasts.fit(zigFit, contrast.matrix)
fit = eBayes(fit)

# Create differential abundance table with every specie identified (datasets 1 and 5 single comparison)
top <-topTable(fit, number = length(tax))
top$OTU <- rownames(top)
top$OTU <- gsub("_", " ", top$OTU)

# Create differential abundance tables for every comparison in datasets 2, 3 and 4 (multiple comparisons, one per coefficient)
topA <-topTable(fit, coef = "a", number = length(tax))
topA$OTU <- rownames(topA)
topA$OTU <- gsub("_", " ", topA$OTU)
topB <-topTable(fit, coef = "b", number = length(tax))
topB$OTU <- rownames(topB)
topB$OTU <- gsub("_", " ", topB$OTU)
topC <-topTable(fit, coef = "c", number = length(tax))
topC$OTU <- rownames(topC)
topC$OTU <- gsub("_", " ", topC$OTU)
topD <-topTable(fit, coef = "d", number = length(tax))
topD$OTU <- rownames(topD)
topD$OTU <- gsub("_", " ", topD$OTU)


# Filter false positives with taxize using taxonomy ranks
# Load libraries
library(taxize)
library(usethis)
usethis::edit_r_environ()                         # Edit environment R file to include NCBI account key needed to connect to the servers
Sys.getenv("ENTREZ_KEY")                          # Check if key needed to connect to NCBI is properly set

# Save all species names from topTable column in a list
b <- top$OTU

# Classificate names to get all taxonomy ranks and tags
taxa_info <- classification(b, db = "ncbi")

# If there are too many taxons (more than 500), taxize tends to fail classification step. 
# It is advised to do this manually in case of unexpected failures during classification.
# Divide list names in shorter lists
b1 <- b[1:500]
b2 <- b[501:1000]
b3 <- b[1001:1500]
b4 <- b[1501:2000]
b5 <- b[2001:2398]

# Classificate every sublist
t1 <- classification(b1, db = "ncbi")
t2 <- classification(b2, db = "ncbi")
t3 <- classification(b3, db = "ncbi")
t4 <- classification(b4, db = "ncbi")
t5 <- classification(b5, db = "ncbi")

# Join lists in a single variable
taxa_info <- c(t1, t2, t3, t4, t5)


# Load species list with every taxonomy filter included (these are exact names that match with the taxonomy names used in NCBI)
# The file includes the tags "Bacteria", "Fungi", "Viruses", "Archaea", species list from EnsemblProtist that will be processed later
filtro <- readLines("/home/alumno/Descargas/TFM_Angel/res2/species2.txt")

# Create empty dataframe to save resulting
df <- data.frame(nombre = character(), taxa = character(), stringsAsFactors = FALSE)

# Loop to process classification results and save matching names between filtro and taxa_info in a new dataFrame
for (i in seq_along(taxa_info)) {
  cla <- taxa_info[[i]]                           # Save iteration in a
  if(any(is.na(cla))){                            # If there is a NA iteration, skip and continue
    next
  }
  else if (any(cla$name %in% filtro)) {           # If there is a name, check if it is in filtro
    noun <- cla$name[cla$rank == "species"]       # Save the name at rank "species" in a variable
    
    if (length(noun) > 0) {                       # If there is something saved, continue, else it will be ignored

      nivel <- filtro[filtro %in% cla$name]       # Asign name used in filtro to the level (Bacteria, Fungi, etc.)
      df <- rbind(df, data.frame(nombre = noun, taxa = nivel))  # Save name and level into a dataframe
    }
  }
  else {                                          # If the iteration is empty, skip and continue
    next
  }
}


# Filter topTable matching OTU names with df taxonomy names
resultado <- subset(top, OTU %in% df$nombre)
# Merge resulting dataFrame with df elements by OTU names to save taxonomy levels (Bacteria, Fungi...)
resultado <- merge(resultado, df, by.x = "OTU", by.y = "nombre", all.x = TRUE)
# Replace all species names from EnsemblProtist with Protist, except the generalist names
resultado$taxa <- ifelse(resultado$taxa %in% c("Bacteria", "Fungi", "Viruses", "Archaea"), resultado$taxa, "Protist")
# Reorder columns for better presentation of tables
resultado <- resultado[, c(2:7, 1, 8)]

# Save topTables in different formats
write.table(resultado[order(resultado$P.Value),],file = "/home/alumno/Descargas/TFM_Angel/res2/final/differential_abundance_filtered/DA_dataset2A_serratia24VSserratia0.tsv", 
          quote = F, sep = "\t", row.names = FALSE, fileEncoding = "UTF-8", dec = ".")
write.table(resultado[order(resultado$P.Value),],file = "/home/alumno/Descargas/TFM_Angel/res2/final/differential_abundance_filtered/DA_dataset2A_serratia24VSserratia0.txt", 
            quote = F, sep = "\t", row.names = FALSE, fileEncoding = "UTF-8", dec = ".")


# Presence-absence tables for comparisons 2B and 2C in dataset 2
# Pre-processing (requires Taxize step to filter)
# Save counts from MRexp normalised object into a variable
t <- MRcounts(mrexp_norm)
colnames(t) <- mdt$Source
t <- t[, c(2, 6, 1, 5, 4, 3)]
t <- as.data.frame(t)
t$OTU <- rownames(t)
t$OTU <- gsub("_", " ", t$OTU) # Run taxize after this step to get the filters

# Filter and save count table (same way than before)
resultado <- subset(t, OTU %in% df$nombre)
resultado <- merge(resultado, df, by.x = "OTU", by.y = "nombre", all.x = TRUE)
resultado$taxa <- ifelse(resultado$taxa %in% c("Bacteria", "Fungi", "Viruses", "Archaea"), resultado$taxa, "Protist")
resultado <- resultado[, c(2:7, 1, 8)]
t <- resultado

# Extract presence-absence information
s0 <- t[, c(1,2)]                            # Serratia 0h counts
rownames(s0) <- t$OTU                        # Replace rownames with OTU names
s24 <- t[, c(3,4)]                           # Serratia 24h counts
rownames(s24) <- t$OTU                       # Replace rownames with OTU names
c0 <- as.data.frame(t[, colnames(t) == "control_0"])     # Control 0h counts
colnames(c0) <- "control_0"                              # Replace colname with control_0
rownames(c0) <- t$OTU                                    # Replace rownames with OTU names
c24 <- as.data.frame(t[, colnames(t) == "control_24"])   # Control 24h counts
colnames(c24) <- "control_24"                            # Replace colname with control_24
rownames(c24) <- t$OTU                                   # Replace rownames with OTU names

# Extract species names that are present in all samples of a condition (2 columns in Serratia and 1 column in Control)
ts0 <- rownames(s0)[rowSums(s0 > 0) > 1]
ts24 <- rownames(s24)[rowSums(s24 > 0) > 1]
tc0 <- rownames(c0)[rowSums(c0 > 0) > 0]
tc24 <- rownames(c24)[rowSums(c24 > 0) > 0]

# Extract taxa that is present in the first condition but absent in the second condition
uniq_s0_c0 <- setdiff(ts0, tc0)             # Taxa present in Serratia 0h that are absent in Control 0h
uniq_c0_s0 <- setdiff(tc0, ts0)             # Taxa present in Control 0h that are absent in Serratia 0h
uniq_s24_c24 <- setdiff(ts24, tc24)         # Taxa present in Serratia 24h that are absent in Control 24h
uniq_c24_s24 <- setdiff(tc24, ts24)         # Taxa present in Control 24h that are absent in Serratia 24h

# Extract dataframe by OTU name with taxa names previously identified
t2 <- t[t$OTU %in% union(uniq_s0_c0, uniq_c0_s0), c(1,2,5,7,8)]         # Taxa from Serratia 0h and Control 0h
t3 <- t[t$OTU %in% union(uniq_s24_c24, uniq_c24_s24), c(3,4,6,7,8)]     # Taxa from Serratia 24h and Control 24h

# Save resulting table in different formats
write.table(t3,file = "/home/alumno/Descargas/TFM_Angel/res2/final/DA_dataset2C_serratia24VScontrol24_presence-absence.txt", 
            quote = F, sep = "\t", row.names = FALSE, fileEncoding = "UTF-8", dec = ".")
write.table(t3,file = "/home/alumno/Descargas/TFM_Angel/res2/final/DA_dataset2C_serratia24VScontrol24_presence-absence.tsv", 
            quote = F, sep = "\t", row.names = FALSE, fileEncoding = "UTF-8", dec = ".")


### PLOTS

# Volcanoplot
# Load libraries
library(ggplot2)
library(ggrepel)
library(scales)

# Load file or use resulting topTable previously created
data <- resultado
data <- read.delim("/home/alumno/Descargas/TFM_Angel/res2/final/DA_dataset2A_serratia24VSserratia0.tsv", header = TRUE, sep = "\t")

# Adapt P-Values to -log values
data$logP <- -log10(data$P.Value)

# Save significantly relevant taxa from literature to represent (run the line of the corresponding comparison)
taxon1 <- c("Plasmodium falciparum", "Exiguobacterium mexicanum", "Nocardia yamanashiensis", "Lachnoclostridium phytofermentans", "Anopheles C virus", "Bacteriophage sp.", "Tsukamurella pulmonis", "Nocardia farcinica", "Serratia marcescens", "Serratia ureilytica", "Paraburkholderia tropica")
taxon2A <- c("Hypocreales sp.", "Saccharomyces cerevisiae", "Serratia marcescens", "Clonostachys rosea", "Asaia bogorensis", "Ophiocordyceps sp.")
taxon2D <- c("Drechmeria coniospora", "Escherichia coli", "Serratia marcescens", "Saccharomyces cerevisiae", "Hypocreales sp.")
taxon3A <- c("Pseudomonas putida", "Anopheles cypovirus", "Acidovorax sp. DW039", "Capnocytophaga ochracea", "Klebsiella pneumoniae", "Bolahun virus variant 1", "Parvimonas micra")
taxon3B <- c("Vibrio vulnificus", "Anopheles C virus", "Actinomyces sp. oral taxon 171", "Bolahun virus variant 1", "Neisseria perflava", "Eimeria tenella")
taxon3C <- c("Neisseria sp. RH3002v2g", "Anopheles C virus", "Sphingomonas sp. FARSPH", "Pseudomonas putida", "Yarrowia lipolytica", "Anopheles cypovirus")
taxon3D <- c("Exophiala lecanii-corni", "Eimeria tenella" ,"Sulfurimonas xiamenesis", "Parvimonas micra", "Candida theae", "Cloacibacterium sp. TD35")
taxon4A <- c("Bradyrhizobium sp.", "Afipia sp.", "Elizabethkingia anophelis", "Bolahun virus variant 1", "Acetobacter aceti", "Pseudomonas tritici", "Aureobasidium melanogenum", "Asaia bogorensis")
taxon4B <- c("Bolahun virus variant 1", "Bradyrhizobium sp.", "Afipia sp.", "Sphingomonas sp. NIBR02145", "Asaia bogorensis", "Delftia lacustris")
taxon4C <- c("Bolahun virus variant 1", "Acetobacter aceti", "Azospirillum baldaniorum", "Aureobasidium melanogenum", "Bradyrhizobium sp. PSBB068", "Penicillium citrinum")
taxon5 <- c("Wallemia sebi", "Bacillus cereus", "Calothrix sp. NIES-4101", "Pseudomonas aeruginosa", "Chryseobacterium piscicola", "Elizabethkingia sp. LA1-18", "Tenacibaculum pacificus")

# Subset the corresponding list in a Top dataframe
Top <- subset(data, OTU %in% taxon3D)
# Shorten tags to present as "E. coli" format. If it is a specie name (sp.) it will not be shortened to keep it legible
Top$OTU <- ifelse(
  grepl(" sp\\.", Top$OTU),
  Top$OTU, sub("^(\\w)\\w*\\s", "\\1. ", Top$OTU)
)

# Prepare tags for easier management of names (content are examples for comparison 3D: Ctrl 48h - Ctrl 24h)
Up <- "Ctrl 48h"
Down <- "Ctrl 24h"
ds <- as.character("3D")
n <- as.character(nrow(data))

# Set variables for proper presentation in volcano
minimoX <- min(data$logFC) - 2
maximoX <- max(data$logFC) + 2
maximoY <- max(data$logP) + 2
down_center <- mean(c(minimoX, -0.58))
up_center <- mean(c(0.58, maximoX))
high <- maximoY - 0.5

# Create volcano plot with differentially abundant taxa
v <- ggplot(data = data, aes(x = logFC, y = logP)) + 
  #Colored areas
  geom_rect(aes(xmin = 0.58, xmax = maximoX, ymin = -log10(0.05), ymax = maximoY), fill = "darksalmon", alpha = 0.01) +
  geom_rect(aes(xmin = minimoX, xmax = -0.58, ymin = -log10(0.05), ymax = maximoY), fill = "cadetblue2", alpha = 0.01) +
  #Points according to filters specified
  geom_point(aes(fill = ifelse(P.Value < 0.05 & logFC > 0.58, "More abundant in Ctrl 48h\n(P-Value < 0.05)", 
                         ifelse(P.Value < 0.05 & logFC < -0.58, "More abundant in Ctrl 24h\n(P-Value < 0.05)", 
                                "Not significant"))), 
             shape = 21, stroke = 0.8, size = 1.5, color = "black") + 
  #Customize color of points previously defined
  scale_fill_manual(values = c("More abundant in Ctrl 48h\n(P-Value < 0.05)" = "red",
                               "More abundant in Ctrl 24h\n(P-Value < 0.05)" = "deepskyblue4", 
                                "Not significant" = "antiquewhite4"), 
                     na.translate = FALSE) + 
  #Title and axis tags
  labs(title = paste0("Differential Abundance Dataset ", ds, ": ", Up, " vs ", Down), x = "Log2 Fold Change", y = "P.Value") + 
  #Legend
  guides(fill = guide_legend(title = paste0("Taxa Identification\n(n = ", n, ")"))) +
  #Simple theme
  theme_minimal() + 
  #Adjust title text positions
  theme(plot.title = element_text(hjust = 0)) +
  #P-Value < 0.05 separating line
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
  
  #Tags for Up, Down and separating line elements
  annotate("text", x = down_center, y = high, label = Down, color = "black", size = 3,
           fontface = "bold") + 
  annotate("text", x = up_center, y = high, label = Up, color = "black", size = 3,
           fontface = "bold") + 
  annotate("text", x = max(data$logFC) + 1, y = -log10(0.05) + 0.5, label = "P = 0.05", hjust = 1, size = 3) +
  
  #Adjust Y axis breaks
  scale_y_continuous(name = "-log P-Value", breaks = seq(0, max(data$logP), by = 5)) +
  #Top taxa
  geom_text_repel(data = Top, aes(label = OTU), size = 3.5, box.padding = 0.5, point.padding = 1, nudge_y = 5, force = 3, 
                  direction = "y", max.overlaps = Inf)
print(v)

# Save in high res svg file
png(paste0("/home/alumno/Descargas/TFM_Angel/res2/final/volcanoplot/complete/Volcano_dataset", ds, "_", Up, "VS", Down, ".png"), 
    width = 4000, height = 4000, res = 600)
print(v)
dev.off()

# Heatmap (presence-absence tables from 2B and 2C)
# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load tables
data <- read.delim("/home/alumno/Descargas/TFM_Angel/res2/final/DA_dataset2B_serratia0VScontrol0_presence-absence.tsv", header = TRUE, sep = "\t")
data <- read.delim("/home/alumno/Descargas/TFM_Angel/res2/final/DA_dataset2C_serratia24VScontrol24_presence-absence.tsv", header = TRUE, sep = "\t")
data1 <- data

# Adjust tables to present a single column with presence or absence. Specie must be in all samples of a condition to be considered present
data1[, -c(4,5)] <- ifelse(data[, -c(4,5)] > 0, 1, 0)
data1[, c(1,2)] <- ifelse(rowSums(data1[, c(1,2)]) < 2, 0, 1)
data1 <- data1[, -2]

# Adapt format to long table to plot heatmap
datalong <- data1 %>% pivot_longer(cols = -c(OTU, taxa), names_to = "sample", values_to = "presence")

# Create heatmap
hm <- ggplot(datalong, aes(x = sample, y = OTU, fill = factor(presence))) + 
  #Use tiles as heat squares
  geom_tile(color = "white") + 
  #Fill tiles with colors according to table values (0 = absence, 1 = present))
  scale_fill_manual(values = c("0" = "grey", "1" = "coral"), 
                    labels = c("0" = "Absent", "1" = "Present"), 
                    name = "Presence") + 
  #Simple theme
  theme_minimal() + 
  #Adjust axis elements to fit and reduce size of heatmap
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = length(unique(datalong$OTU)) / length(unique(datalong$sample))) + 
  #Title and axis tags
  labs(title = "Heatmap Presence in Serratia 24h vs Control 24h", x = "Groups", y = "Species")
print(hm)

# Save heatmap in a high res svg file
png("/home/alumno/Descargas/TFM_Angel/res2/final/Heatmap_dataset2C_serratia24VScontrol24_presence-absence.png", width = 5000, height = 6000, res = 600)
print(hm)
dev.off()