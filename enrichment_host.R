# Parte Enrichment (edgeR)
library("edgeR")
# Lee tabla y crea factores
setwd("/home/alumno/miniconda3/envs/qiime2-dev/kraken2")
data <- read.delim("fuse_matrixes/htseq-count_PE.tsv", header = TRUE, row.names = 1)
data <- read.delim("fuse_matrixes/htseq-count_subsampled_PE.tsv", header = TRUE, row.names = 1)
data <- read.delim("fuse_matrixes/htseq-count_new_subsampled_PE.tsv", header = TRUE, row.names = 1)
head(data)
factores <- data.frame(Source = c("midgut", "midgut", "salivary-gland", "salivary-gland"), row.names = colnames(data))
factores$Source <- factor(factores$Source, levels = c("salivary-gland", "midgut"))
#group <- c("midgut", "midgut", "salivary-gland", "salivary-gland")
# Factores para muestras nuevas
factores <- data.frame(Type = c("infected", "infected", "control", "control"), row.names = colnames(data))
factores$Type <- factor(factores$Type, levels = c("infected", "control"))


# Crea objeto DGE, ordena y comprueba duplicados
dge <- DGEList(counts = data, group = factores$Source, genes = rownames(data))
# DGE para muestras nuevas
dge <- DGEList(counts = data, group = factores$Type, genes = rownames(data))
ord <- order(rowSums(dge$counts), decreasing = TRUE)
dge <- dge[ord,]
dupes <- which(duplicated(dge$genes))
length(which(dupes == TRUE))

# Comprueba balanceo
dge$samples$lib.size <- colSums(dge$counts)
(dge$samples$lib.size - mean(dge$samples$lib.size)) / mean(dge$samples$lib.size) * 100  # Muestra 4 NO ESTA BALANCEADA

# Downsampling de muestras 1, 2 y 3 segun formula. Ej: 
# 1) 14335359 -> 21517704 reads, luego para 10892209, necesitaria X
10892209*21517704/14335359
# 2) 14215520 -> 25932154 reads
10892209*25932154/14215520
# 3) 12809434 -> 18753987 reads
10892209*18753987/12809434

# Modelo (design)
modelo <- model.matrix(~factores$Source)
#modelo <- model.matrix(~factores$Source + factores$Number)
rownames(modelo) <- rownames(factores)

# Modelo para muestras nuevas
modelo <- model.matrix(~factores$Type)
rownames(modelo) <- rownames(factores)

# Filtro de expresion por defecto y normalizado
keep <- filterByExpr(dge, modelo)                               # Filtro expresión retiene 8074 genes de 13850 / 7785 subsampled
table(keep)
dge <- dge[keep,]
dge <- normLibSizes(dge)                                        # Normalizado
#filterByExpr retiene solo 629 genes de los 13850

# Estima dispersion y crea modelo final
dge <- estimateDisp(dge, modelo, robust = TRUE)
fit <- glmFit(dge, modelo)
lrt <- glmLRT(fit)
topTags(lrt)                                                    # Comprueba topTags (no se usa, orientativo)

# Ajuste valor p (FDR)
lrt$table$FDR <- p.adjust(lrt$table$PValue, method = "BH")
# Filtra por FDR < 0.05
genes <- rownames(lrt$table)
genesFiltered <- lrt$genes[lrt$table$FDR < 0.05,]
length(genes)
length(genesFiltered)                                                   # Suma total de genes up y down = 2868

# Ordena por FDR (arreglar)
ord <- order(lrt$table$FDR)
lrt$table <- lrt$table[ord,]

# Decide genes y volcano plot apañao
#cpm(dge)[ord[1:10],]
summary(decideTests(lrt))                                                 # Test para ver genes up y down
plotMD(lrt)
abline(h=c(-1, 1), col="green")                                           # Volcano plot apañado

# Guardar tabla 2868 genes (FDR < 0.05)
x <- lrt$table[which(lrt$table$FDR < 0.05 & lrt$table$logFC < 0),]
x <- data.frame(VB_ID = rownames(x), x)
write.table(rownames(x), file = "genelist_down.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(x, file = "deg_downsampling.tsv", sep = "\t", row.names = FALSE)

#Coger suma genes up-down y hacer enrichment con universe = todos los genes detectados
x <- lrt$table[which(lrt$table$FDR < 0.05),]
length(genes)
length(genesFiltered)
lfc <- 1
padj <- 0.05
up <- x[which(x$logFC > lfc & x$FDR < padj), ]
up <- rownames(up)
lfc <- -1
down <- x[which(x$logFC < lfc & x$FDR < padj), ]
down <- rownames(down)
length(up)
length(down)


# Enricher (incompleta)
library(clusterProfiler)
gaf <- read.delim("gff files/VectorBase-67_AgambiaePEST_Curated_GO.gaf", header = FALSE, comment.char = "!", stringsAsFactors = FALSE)
gene2GO <- gaf[, c(2, 5)]
colnames(gene2GO) <- c("gene", "GO")
go_notes <- split(gene2GO$GO, gene2GO$gene)
backgroundGO <- unique(gene2GO$gene)
background <- which(genes %in% backgroundGO)
backgroundF <- genes[background]

enrichmentUp <- enricher(gene = up, universe = backgroundF, TERM2GENE = gene2GO[,c(2, 1)])
enrichmentDown <- enricher(gene = down, universe = backgroundF, TERM2GENE = gene2GO[,c(2, 1)])
resUp <- as.data.frame(enrichmentUp)
resDown <- as.data.frame(enrichmentDown)

#Coger GOs y darles nombre
library(GO.db)
resUp$Description <- Term(resUp$ID)
resDown$Description <- Term(resDown$ID)

#GSEA
ord <- order(lrt$table$logFC, decreasing = TRUE)
lord <- lrt$table[ord,]

listado <- lord$logFC
names(listado) <- rownames(lord)
listado <- sort(listado, decreasing = TRUE)

gseaenrich <- GSEA(geneList = listado, TERM2GENE = gene2GO[, c(2, 1)], pvalueCutoff = 0.05)
resGSEA <- as.data.frame(gseaenrich)

library(GO.db)
resGSEA$Description <- Term(resGSEA$ID)

########
y <- x$logFC
names(y) <- rownames(x)
# Añadir TERM2NAME
enrico <- enricher(gene = rownames(x), TERM2GENE = data.frame(GO = gene2GO$GO, gene = gene2GO$gene), 
                   pvalueCutoff = 0.05, universe = rownames(lrt$table))
res <- as.data.frame(enrico)[,c(2:7,9)]
Term <- c("serine-type endopeptidase activity ", "extracellular space", "proteolysis", "extracellular region", "plasma membrane", "serine-type peptidase activity", "transmembrane transport", "peptidase activity", "plasma membrane")
write.table(res, file = "table_enricher.tsv", sep = "\t", row.names = FALSE)

gsearico <- GSEA(geneList = sort(y, decreasing = TRUE), TERM2GENE = data.frame(GO = gene2GO$GO, gene = gene2GO$gene),
                pvalueCutoff = 0.05, universe = rownames(lrt$table))

#Notas
# Nunca dividir ups y downs al hacer enrichment
# Comprobar si muestras estan balanceadas (conteo +-10% maximo -> downsampling)
# Cuando pida universo, coger todos genes detectados en experimentos
# Usar filtro expresion (filterByExpr) como universo
# Normalizar y DEG
# Usar AnnotationForge para construir OrgDb, pero se supone que el Org.ag de anopheles esta actualizado -> No funciona


# Alternativa con topGO (funciona, la mas sencilla)
library("topGO")
gaf <- read.delim("gff files/VectorBase-67_AgambiaePEST_Curated_GO.gaf", header = FALSE, comment.char = "!", stringsAsFactors = FALSE)
gene2GO <- gaf[, c(2, 5)]
colnames(gene2GO) <- c("gene", "GO")
geneID2GO <- split(gene2GO$GO, gene2GO$gene)
#backgroundGO <- unique(gene2GO$gene)                     # Por si hay dupes

#Identifica GOs con los genes deg
#allGenes
geneList <- factor(as.integer(genes %in% genesFiltered))
names(geneList) <- genes

#upGenes
upGenes <- factor(as.integer(genes %in% up))
names(upGenes) <- genes

#downGenes
downGenes <- factor(as.integer(genes %in% down))
names(downGenes) <- genes

#ggplot individual (GO compara allGenes como factor de si está o no del total de background)
ont <- "BP"
GOdata <- new("topGOdata", 
              ontology = ont, 
              allGenes = upGenes, 
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

# Extrae tabla enrichment y valores importantes
enrichment <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata, classicFisher = enrichment, topNodes = 50, numChar = 100)
allRes$classicFisher <- as.numeric(gsub("<", "", allRes$classicFisher))               # Si el valor es demasiado pequeño, elimina "<"
allRes$adjFisher <- p.adjust(allRes$classicFisher, method = "BH")
#allRes$geneRatio <- allRes$Annotated / allRes$Significant
#allRes <- allRes[which(allRes$Fisher.adjust < 0.01), ]
#allRes$Annotated <- as.numeric(allRes$Annotated)
#allRes$Significant <- as.numeric(allRes$Significant)
#allRes$Expected <- as.numeric(allRes$Expected)
#allRes$classicFisher <- as.numeric(allRes$classicFisher)
#allRes$Fisher.adjust <- as.numeric(allRes$Fisher.adjust)
print(allRes)


###################
a <- topTags(lrt, n = 100)
logfc_threshold <- 0.5
padj_threshold <- 0.05
up <- x$table$genes[which(x$table$logFC > logfc_threshold & x$table$PValue < padj_threshold)]
down <- x$table$genes[which(x$table$logFC < logfc_threshold & x$table$PValue < padj_threshold)]
x <- lrt$genes[which(lrt$table$FDR < padj_threshold),]      # Test

# Bucle para genes UP y DOWN y sacar ggplot de todas las ONT
for (n in c("Down", "Up")) {
  #deg <- down
  pdf(paste("Enrichment_", n, "_classicFisher.pdf", sep = ""), paper = "A4")
  logfc_threshold <- 0
  padj_threshold <- 0.05
  
  if(n == "Down"){
    deg <- lrt$table[which(lrt$table$logFC < logfc_threshold & lrt$table$PValue < padj_threshold)]
  }else{
    deg <- x$table$genes[which(x$table$logFC > logfc_threshold & x$table$PValue < padj_threshold)]
  }
  
  geneList <- factor(as.integer(background %in% deg))          # Puede ir DEgenes, up or down genelists
  names(geneList) <- background
  
  for (ONT in c("BP", "MF", "CC")) {

    GOdata <- new("topGOdata", 
                  ontology = ONT, 
                  allGenes = geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = geneID2GO)
    enrichment <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    allRes <- GenTable(GOdata, classicFisher = enrichment, topNodes = 50, numChar = 1000)
    Fisher.adjust <- p.adjust(allRes$classicFisher, method = "BH")
    allRes <- cbind(allRes, Fisher.adjust)
    print(allRes)
    # Hasta aquí, lista enrichment que FUNCIONA
    
    # Gráfica enrichment
    library("ggplot2")
    library("RColorBrewer")
    
    cond <- as.numeric(allRes$Fisher.adjust) < 0.05
    tt <- allRes[cond,]
    tt$Annotated <- as.numeric(tt$Annotated)
    tt$Significant <- as.numeric(tt$Significant)
    tt$Expected <- as.numeric(tt$Expected)
    tt$classicFisher <- as.numeric(tt$classicFisher)
    tt$Fisher.adjust <- as.numeric(tt$Fisher.adjust)
    if (length(tt$GO.ID) < 5){
      tt <- allRes[1:5,]
    }
    paleta <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
    ordenes <- 10^-(0:30)
    
    g <- ggplot(data = tt, aes(x = reorder(Term, Fisher.adjust, decreasing = TRUE), y = Significant)) +
      geom_bar(stat = "identity", color = "black", aes(fill = log(Fisher.adjust)), size = 0.3) + 
      geom_text(aes(label = Term), color = "black", position = position_fill(vjust = 0), hjust = 0, 
                fontface = "bold", size = 4) +
      coord_flip() + 
      theme(
        panel.background = element_blank(), 
        panel.grid.major.x = element_line(colour = "darkgrey", linewidth = 0.75), 
        panel.grid.minor.x = element_line(colour = "grey", linewidth = 0.75), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.y = element_blank(), 
        axis.text = element_text(size = 12)) + 
      ylab(paste("Significant Genes", n, ONT)) + 
      guides(fill = guide_colourbar(barheight = 25, reverse = T)) + 
      scale_fill_gradientn(name = "Fisher-adjusted", 
                           colours = paleta(4), 
                           breaks = log(ordenes), 
                           guide = guide_colourbar(reverse = TRUE), 
                           labels = ordenes) +
      scale_y_continuous(breaks = seq(0, max(tt$Significant), by = 5))
    print(g)
  }
  dev.off()
}

library("ggplot2")
library("RColorBrewer")

paleta <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
ordenes <- 10^-(0:30)

g <- ggplot(data = allRes, aes(x = reorder(Term, Fisher.adjust, decreasing = TRUE), y = geneRatio)) +
  geom_bar(stat = "identity", color = "black", aes(fill = as.numeric(Fisher.adjust)), size = 0.3) + 
  geom_text(aes(label = Term), color = "black", position = position_fill(vjust = 0), hjust = 0, 
            fontface = "bold", size = 4) +
  coord_flip() + 
  theme(
  panel.background = element_blank(), 
  panel.grid.major.x = element_line(colour = "darkgrey", linewidth = 0.75), 
  panel.grid.minor.x = element_line(colour = "grey", linewidth = 0.75), 
  axis.title.y = element_blank(), 
  axis.text.y = element_blank(), 
  axis.ticks.y = element_blank(), 
  axis.ticks.x = element_blank(), 
  axis.line.y = element_blank(), 
  axis.text = element_text(size = 12)) + 
  ylab(paste("Significant Genes")) + 
  guides(fill = guide_colourbar(barheight = 25, reverse = T)) + 
  scale_fill_gradientn(name = "Fisher.adjust", 
                       colours = paleta(9), 
                       breaks = log(ordenes), 
                       guide = guide_colourbar(reverse = TRUE), 
                       labels = ordenes) +
  scale_y_continuous(breaks = seq(0, max(allRes$geneRatio)))
print(g)

## Añadir ajuste Fisher tipo BH para tabla resultado
weight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
allRes1 <- GenTable(GOdata, weightFisher = weight, orderBy = "weightFisher", topNodes = 50)
p.adj1 <- round(p.adjust(allRes$classicFisher, method = "BH"), digits = 4)
allRes2 <- cbind(allRes1, p.adj1)


# Otras graficas
plotMDS(dge, col=rep(1:2, each=2))

# Dispersion
dge$common.dispersion
plotBCV(dge)

# Otra dispersion
fit2 <- glmQLFit(dge, modelo, robust = TRUE)
plotQLDisp(fit2)



# Parte Enrichment (DESeq2, falla Org.Db)
library("DESeq2")
setwd("/home/alumno/miniconda3/envs/qiime2-dev/kraken2")
data <- read.delim("fuse_matrixes/htseq-count_PE.tsv", header = TRUE, row.names = 1)
head(data)
factores <- data.frame(Source = c("midgut", "midgut", "salivary", "salivary"), 
                       Number = c("1", "2", "1", "2"), row.names = colnames(data))
factores$Source <- factor(factores$Source, levels = c("midgut", "salivary"))

dds <- DESeqDataSetFromMatrix(countData = as.matrix(data),
                              colData = DataFrame(factores), 
                              design = ~ Source + Number)

model.matrix(~Source + Number, data = factores)

keep <- rowSums(counts(dds)) >= 10  # No funciona? No elimina ninguno?
dds2 <- dds[keep,]

res <- DESeq(dds2, test = "Wald")   # Retiene 9654 genes
DESeq2::plotSparsity(res)
DESeq2::plotDispEsts(res)

vst_dds <- vst(res)
DESeq2::plotPCA(vst_dds, intgroup = c("Source"), ntop = nrow(res))

resultsNames(res)
res_trt <- results(res, name = "Source_salivary_vs_midgut", pAdjustMethod = "fdr")
logfc_threshold <- 0
padj_threshold <- 0.05

table(res_trt$log2FoldChange > logfc_threshold & res_trt$padj < padj_threshold) # Genes sobre-expresados 1572
table(res_trt$log2FoldChange < logfc_threshold & res_trt$padj < padj_threshold) # Genes sub-expresados   867


# Anotacion
library("AnnotationDbi")
library("org.Ag.eg.db")

columns(org.Ag.eg.db)
rownames(res_trt)
head(res_trt)

res_trt$ID <- rownames(res_trt)

res_trt$symbol <- mapIds(org.Ag.eg.db, keys = rownames(res_trt), column = "SYMBOL", keytype = "ENSEMBL")

res_trt$GO <- NA
for (i in 1:nrow(gaf)) {
  gene <- res_trt$ID[i]
  go_terms <- gaf$GO[gaf$GENE == gene]
  res_trt$GO[i] <- paste(go_terms, collapse = ", ")
}

gaf2 <- gaf[, c("GENE", "GO")]
gaf2 <- unique(gaf2[,1])
res_trt2 <- merge(res_trt, gaf2, by = "GENE", all.x = TRUE)

gaf <- read.gaf("gff files/VectorBase-67_AgambiaePEST_Curated_GO.gaf")
gaf <- gaf[order(gaf$GENE),]
agap_genes <- subset(gaf, grepl("^AGAP", gaf$GENE))
go_terms <- agap_genes$GO


logfc_threshold <- 0.5
padj_threshold <- 0.05
deg <- list(up = res_trt$ID[which(res_trt$log2FoldChange > logfc_threshold & res_trt$padj < padj_threshold)], 
            down = res_trt$ID[which(res_trt$log2FoldChange < -logfc_threshold & res_trt$padj < padj_threshold)])

# Hacer clusterProfiller para enrichment (o coger topGO)

library("clusterProfiler")

genetoGO <- split(gaf$GO, gaf$GENE)
x <- enricher(gene = deg$up, TERM2GENE = , TERM2NAME = )

fea_GO_BP <- compareCluster(deg, ont = "BP", fun = "enrichGO", idx = genetoGO)
fea_GO_MF <- compareCluster(deg, ont = "MF", fun = "enrichGO")
fea_GO_CC <- compareCluster(deg, ont = "CC", fun = "enrichGO")


