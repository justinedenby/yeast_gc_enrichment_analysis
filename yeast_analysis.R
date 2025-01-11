library(RERconverge)

#yeast master tree
yeastTree <- readRDS("yeast_tree.tre.RDS")
names(yeastTree$trees) <- c(new_vector) #this renames 
summary(yeastTree)

plot(yeastTree$masterTree)

# RER function
yeastRER = getAllResiduals(yeastTree, transform = "sqrt", weighted = T, scale = T)
saveRDS(yeastRER, file="yeastRER.rds")

#two data sets for continuous trait analysis
yeast_RER = readRDS("yeastRER.rds")
gc_contentdf = read.csv("Yeast_GC_Percents.csv")
gc_content<- gc_contentdf$GC.Content
names(gc_content) <- gc_contentdf$File
print(gc_content)

# gene names to a variable
geneNames = read.table("FullScerGeneOGConversion.tsv", header=TRUE)

#gene renaming
OGvector <- (names(yeastTree$trees))
names(OGvector) <- OGvector
GeneVector <- geneNames$YCF1
names(GeneVector) <- geneNames$OG1000

# finds matching values and replace names with values from GeneVector
new_vector <- OGvector
matching_indices <- match(OGvector, geneNames$OG1000) 
matching_values <- !is.na(matching_indices)  # checks which values have a match
new_vector[matching_values] <- geneNames$YCF1[matching_indices[matching_values]]

# continuous trait analysis, RER specific
charP = char2Paths(gc_content, yeastTree)
hist(charP)

#correlate phenotype (gc content), RER specific
cphen = correlateWithContinuousPhenotype(yeast_RER, charP, min.sp = 10,
                                         winsorizeRER = 2, winsorizetrait = 0)
hist(cphen$P, breaks=35)
#Rho, p-value, and Benjamini-Hochberg corrected p-value

cphen$stat = -log10(cphen$P) * sign(cphen$Rho)
head(cphen[order(cphen$stat, decreasing = TRUE),])
#these are top positively correlated genes. if you change decreasing
#to false, it lists the negative correlations

####################################

#pathway gene enrichment
library(RERconverge)
data("RERresults")
print(data)
stats=getStat(RERresults)
# RER specific

annots1=read.gmt("GO_Molecular_Function_2018.txt")
annots2=read.gmt("WikiPathways_2018.txt")
annots3=read.gmt("KEGG_2019.txt")
annots4=read.gmt("GO_Cellular_Component_2018.txt")
annots5=read.gmt("GO_Biological_Process_2018.txt")
annots6=read.gmt("Phenotype_AutoRIF.txt")

vals = getStat(cphen)
#builds object to hold enrichment data
annotlist=list(annots1,annots2,annots3,annots4,annots5,annots6)
names(annotlist) <- c("GOMolecular","WikiPathways","KEGG","GOCellular","GOBiological", "Phenotype")

enrichment=fastwilcoxGMTall(vals, annotlist, outputGeneVals=T, num.g=10)
View(enrichment)

####

# view the stat, pval, and p.adj of the top enrichment results
head(enrichment$GOMolecular[order(enrichment$GOMolecular$pval),])[1:3]
View(enrichment$GOMolecular)
print(enrichment$GOMolecular[order(enrichment$GOMolecular$pval),])[1:3]

#Wiki Pathways
head(enrichment$WikiPathways[order(enrichment$WikiPathways$pval),])[1:3]
View(enrichment$WikiPathways)

#KEGG
head(enrichment$KEGG[order(enrichment$KEGG$pval),])[1:3]
View(enrichment$KEGG)

#GO Cellular
head(enrichment$GOCellular[order(enrichment$GOCellular$pval),])[1:3]
View(enrichment$GOCellular)

#GO Biological
head(enrichment$GOBiological[order(enrichment$GOBiological$pval),])[1:3]
View(enrichment$GOBiological)

#Phenotype
head(enrichment$Phenotype[order(enrichment$Phenotype$pval),])[1:3]
View(enrichment$Phenotype)