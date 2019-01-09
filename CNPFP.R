# read in the R libraries
library(MASS)	# standard, no need to install
library(class)	# standard, no need to install
library(cluster)	# standard, no need to install
library(sma)	# install it for the function plot.mat 
library(impute)	# install it for imputing missing value
library(Hmisc)	# install it for the C-index calculations
# read in the custom network functions
source("(file location)/NetworkFunctions.txt")
library(WGCNA)
options(stringsAsFactors = F)

##Network Construction and ANALYSIS OF YPPI & YGE DATASETS
source("(file location)/YPPINetworkFunctions.txt");
source("(file location)/YPPINetworkFunctions1.txt");
NodeNames=read.csv(file = "(file location)/YeastProtein_NodeNames.csv", header=F)
NodeNames=as.vector(NodeNames[,1])
ADJ1=read.csv("(file location)/YeastProtein_ADJ1.csv", header=F)
ADJ1=as.matrix(ADJ1)
dim(ADJ1)
rownames(ADJ1)=NodeNames; 
colnames(ADJ1)=NodeNames;
datExprADJ1 <- as.data.frame(t(ADJ1))
names(datExprADJ1)
# We use the TOM matrix as a smoothed-out version of the adjacency matrix.
dissGTOM1=TOMdist1(as.matrix(ADJ1))
ADJ=1- dissGTOM1
Degree <- apply(ADJ,2,sum)
DegCut = 2292 # number of most connected genes that will be considered 
DegreeRank <- rank(-Degree)	
restDegreeADJ <- DegreeRank <= DegCut
sum(restDegreeADJ)
AdjMat1rest <- ADJ[restDegreeADJ,restDegreeADJ]
dissGTOM1=TOMdist1(as.matrix(AdjMat1rest))
collect_garbage()

# Now we carry out hierarchical clustering with the TOM matrix. Branches of the 
# resulting clustering tree will be used to define gene modules.

hierGTOM1 <- hclust(as.dist(dissGTOM1),method="average");
windows(width=8, height=5)
par(mfrow=c(1,1), mar=c(0,2,2,0))
plot(hierGTOM1,labels=F, cex=0.2)
myheightcutoff    =0.97
mydeepSplit       =  T # fine structure within module
myminModuleSize   = 50 # modules must have this minimum number of genes
myminAttachModuleSize=5
#new way  for identifying modules based on hierarchical clustering dendrogram
colorh1=cutreeDynamic(hierclust= hierGTOM1, deepSplit=mydeepSplit, maxTreeHeight=myheightcutoff, minModuleSize=myminModuleSize, minAttachModuleSize=myminAttachModuleSize, nameallmodules=T, useblackwhite=F)

table(colorh1)
colorlevel1=levels(factor(colorh1))

#Hierarchical clustering dendrogram and module definition
windows()
par(mfrow=c(2,1), mar=c(2,2,2,1))
plot(hierGTOM1, main="Yeast Protein-protein Interaction Network", labels=F, xlab="", sub="");
hclustplot1(hierGTOM1,colorh1, title1="",colorHeight = 0.1)

collect_garbage()
pdf("DendrogramOfGenes.pdf", width=15, height=8)
plotDendroAndColors(hierGTOM1,colorh1,"",main= "Dendrogram of modules", dendroLabels=FALSE)
dev.off();

windows()
par(mfrow=c(1,2))
plotDendroAndColors(hierGTOM1, colorh1, "", autoColorHeight = FALSE,
                    colorHeight = 0.2,
                    ylab = "", xlab = "", sub = "", dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, cex.rowText = 1.3, cex.colorLabels = 1.2,
                    guideHang = 0.05, 
                    main = "Yeast Protein-protein Interaction module colors")

###############Detection of hub genes
ConnectivityMeasuresdataOne=intramodularConnectivity(AdjMat1rest,colors=colorh1) 
names(ConnectivityMeasuresdataOne) 
colorlevels=levels(factor(colorh1)) 
par(mfrow=c(2,3)) 
for (i in c(1:length(colorlevels) ) ) { 
  whichmodule=colorlevels[[i]];restrict1=colorh1==whichmodule 
  verboseScatterplot(ConnectivityMeasuresdataOne$kWithin[restDegreeADJ], AdjMat1rest,col=colorh1,main= paste("set I,", whichmodule),ylab="Gene Significance",xlab="Intramodular k") 
} 

# The following verifies that module definition is relatively robust to the height cutoff.
# This is especially true for tighter modules (more distinctive branches).
windows(width=8, height=10)
par(mfrow=c(6,1), mar=c(2,2,2,1))
plot(hierGTOM1, main="Standard TOM Measure", labels=F, xlab="", sub="");
hclustplot1(hierGTOM1,colorh1, title1="Our chosen height cutoff = 0.93")
for(varheightcutoff in c(0.92, 0.94, 0.95, 0.96) ){
  colorhvar=cutreeDynamic(hierclust= hierGTOM1, deepSplit=mydeepSplit, maxTreeHeight=varheightcutoff, minModuleSize=myminModuleSize, minAttachModuleSize=myminAttachModuleSize, nameallmodules=T, useblackwhite=F)
  hclustplot1(hierGTOM1,colorhvar, title1=paste("Height cutoff =", varheightcutoff)) 
}

#Calculate the fundamental network concepts
colorlevel1=levels(factor(colorh1))
grey=45 # To flag which module is the grey module
ADJ=as.matrix(ADJ)

#Heatmap for all modules
# for visualizing the network. Here we chose 2 scaling dimensions
windows()
cmd1=cmdscale(as.dist(dissGTOM1),2)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorh1),  main="MDS plot", xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

windows()
TOMplot1(dissGTOM1 , hierGTOM1, colorh1)

# Calculate eigengenes 
require(WGCNA)
PMEList = moduleEigengenes(datExprADJ1[restDegreeADJ], colors = colorh1) 
PMEs = PMEList$eigengenes 
# Calculate dissimilarity of module eigengenes 
PMEDiss = 1-cor(PMEs); 
# Cluster module eigengenes 
PMETree = hclust(as.dist(PMEDiss), method = "average"); 
# Plot the result 
windows() 
plot(PMETree, main = "Clustering of module eigengenes", xlab = "", sub = "")


##Module membership
geneModuleMembership1 = signedKME(datExprADJ1[restDegreeADJ], PMEs)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep="")
MMPvalue1 = corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExprADJ1[restDegreeADJ])[[2]]);
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep="");
Gene = rownames(datExprADJ1[restDegreeADJ])
kMEtable1 = cbind(Gene,Gene,colorh1)
for (i in 1:length(colorsA1))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i],
                    MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(
  geneModuleMembership1), colnames(MMPvalue1))))

#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge 
MEDissThres = 0.25 
# Call an automatic merging function 
merge = mergeCloseModules(dissGTOM1, colorh1, cutHeight = MEDissThres, verbose = 3) 
# The merged module colors 
mergedColors = merge$colors;
# Eigengenes of the new merged modules: 
mergedMEs = merge$newMEs;
#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
windows()
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6) 
plotDendroAndColors(hierGTOM1, cbind(colorh1, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
# Rename to moduleColors 
moduleColors = mergedColors 
# Construct numerical labels corresponding to the colors 
colorOrder = c("grey", standardColors(50)); 
moduleLabels = match(moduleColors, colorOrder)-1;
table(moduleLabels)
moduleColors = labels2colors(colorh1)
MEs = mergedMEs; 
geneTree = merge$dendrograms[[1]];


###Yeast Cell Cycle dataset analysis
library(WGCNA)
source("(file location)/YGENetworkFunctions.txt"); 
source("(file location)/YGENetworkFunctions1.txt");
# Prepare the data
dat0 <- read.csv(file ="(file location)/YEASTCellCycle4000.csv",header=T, row.names=1)
dim(dat0)
datExpr = t(dat0[,8:51])
names(datExpr)
datExprYeast=as.data.frame(t(dat0[,-c(1:7)]))
names(datExprYeast)
Yeastdata <- na.omit(dat0)
dim(Yeastdata)
datExprYeast1=as.data.frame(t(Yeastdata[,-c(1:7)]))
str(datExprYeast1)

# This following code allow us to restrict our analysis
# to the most varying genes
var1=function(x) var(x,na.rm=T)
vargenes=apply(datExpr,2,var1)
rankvar=rank(-vargenes)
restVariance= rankvar< 1264 #i.e. we keep all genes.
sum(restVariance)
datExpr=datExpr[,restVariance]
datSummary=datSummary[restVariance,]
no.samples = dim(datExpr)[[1]]
dim(datExpr)
rm(Yeastdata);
gc()


require(graphics)
require(flashClust)
sampleTree = flashClust(dist(datExpr), method = "average");
windows()
par(cex = 0.9);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex = 1.5,cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# soft thresholds (powers). It will be used to pick a soft threshold.
powers = c(seq(from = 1, to=20, by=1))
powers1=c(seq(1,10,by=1), seq(12,20, by=2))
# Call the network topology analysis function
#Using maxPOutliers=1 will effectively disable all weight function broadening
require(WGCNA)
sft = pickSoftThreshold(datExprYeast1, powerVector = powers, verbose = 5, networkType ="signed hybrid", corFnc= "bicor", corOptions=list(maxPOutliers=0.1))
# Scale-free topology fit index as a function of the soft-thresholding power
library("ggplot2")
library("gplots")
library("RColorBrewer")
windows()
par(cex = 0.6);
par(mfrow=c(1,2))
par(mar = c(0,4,2,0))
ggplot(data = sft$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = 0.9), colour = "red")

ggplot(data = sft$fitIndices, aes(x = Power, y = mean.k.)) +
  geom_point(size = 3) +
  geom_line(size = 0.5) +
  ggtitle(label = "Mean connectivity")

# Mean connectivity as a function of the soft-thresholding power
softPower = 5
cor.Y <- bicor(datExprYeast1, use = "pairwise.complete.obs", maxPOutliers = 0.1)
adj = adjacency.fromSimilarity(cor.Y, type = "signed hybrid", power = softPower)
TOM = TOMsimilarity(adj, TOMDenom = "mean", TOMType = "signed")
colnames(TOM) <- rownames(TOM) <- rownames(Yeastdata)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM),method="average");
#plot the resulting clustering tree (dendrogram)
windows()
cex1=1
par(mfrow=c(1,2))
plot(geneTree, xlab="", sub="",cex=0.1, ann = F);
myheightcutoff    =0.99
mydeepSplit       =  TRUE # fine structure within module
myminModuleSize   = 50 # modules must have this minimum number of genes
# Set the minimum module size
#minModuleSize = 30;
# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(hierclust = geneTree,deepSplit=mydeepSplit,maxTreeHeight=myheightcutoff,minModuleSize=myminModuleSize)
table(dynamicMods) 


###Hub Genes 
hubs = chooseTopHubInEachModule(datExprYeast1, dynamicMods)
hubs
  
  ## Make an adjacency matrix and set some arbitrary cutoffs and parameters for plotting the network
  c = as.vector(data)
  y <- subset(datExprYeast1, rownames(datExprYeast1)%in%c)
  adx <- as.data.frame(t(y[,-c(3:51)]))
  addj <- as.matrix(y)
  addj <- bicor(y, use = "pairwise.complete.obs", maxPOutliers = 0.1)
  #adjMat <- bicor(datExprYeast1, use = "pairwise.complete.obs", maxPOutliers = 0.1)
  softPower = 5
  cor.Y <- bicor(datExprYeast1, use = "pairwise.complete.obs", maxPOutliers = 0.1)
  adjMat = adjacency.fromSimilarity(cor.Y, type = "signed hybrid", power = softPower)
  maxsize <- min(70,nrow(adjMat))
  metsize <- min(30,nrow(adjMat))
  adjMat[adjMat<0.2] <- 0
  kME <- apply(adjMat,2,sum)
  cutoff <- sort(kME,decreasing=TRUE)[maxsize]
  cutoff2 <- sort(kME,decreasing=TRUE)[metsize]
  tophubs <- kME>=cutoff
  metgenes <- names(kME)[kME>=cutoff2]
  keepgenes <- NULL
  keepgenes <- c(keepgenes,metgenes)
  adjMat <- adjMat[tophubs,tophubs]
  numcors <- min(200,(maxsize^2-maxsize)/2)
  topcors <- sort(as.numeric(adjMat),decreasing=TRUE)[numcors]
  adjMat[adjMat<=topcors] <- 0
  
  ## Create a network object to store the data, plot these networks in a circle plot
  gA <- graph.adjacency(as.matrix(adjMat[1:5,1:5]),mode="undirected",weighted=TRUE,diag=FALSE) ## Top 5 in the center
  gB <- graph.adjacency(as.matrix(adjMat[6:maxsize,6:maxsize]),mode="undirected",weighted=TRUE,diag=FALSE) ## Additioanl genes on the periphery
  layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB)) ## Construct layout
  
  g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE) ## Graph information for plotting
  
  tkid <- tkplot(g1)
  l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
  
  tk_close(tkid, window.close = T)
  windows()
  plot(g1,vertex.label=geneSymbols,
       vertex.label.dist=0.3, 
       vertex.shape="none",
       vertex.size=6,
       vertex.label.color="black", 
       vertex.label.cex=.7, 
       vertex.color="orange",
       vertex.label.cex=0.8, vertex.label.dist=3,
       vertex.frame.color="black",
       layout=l,
       edge.color="Black")
  
  graph.density(g1)
  
  cor(authority.score(g1)$vector, V(g1)$auth)
  cor(hub.score(g1)$vector, V(g1)$hub)
  hs <- hub_score(g1, weights=NA)$vector
  as <- authority_score(g1, weights=NA)$vector
  Degree <- degree(g1)
  centralities <- cbind(Degree, hs, as)
  write.csv(centralities, file="centralitiesPPI.csv")
  cor(centralities)
  Betweenness <- betweenness(g1)
  blocks <- cohesive.blocks(g1)
  blocks
  summary(blocks)
  par(mfrow=c(1,2))
  windows()
  plot(g1, vertex.size=hs*50, main="Hubs")
  windows()
  plot(g1, vertex.size=as*30, main="Authorities")

library(WGCNA)
library(igraph)
pow=6
net.1 = blockwiseModules(datExpr.1, power = pow,
                         maxBlockSize = 10000, deepSplit = 2,
                         minModuleSize = 10,
                         saveTOMs = FALSE,
                         verbose = F)
module_assignments <- net.1$colors
nlevels(factor(module_assignments))

hs <- hub_score(dynamicMods, weights=NA)$vector

as <- authority_score(net, weights=NA)$vector
## Here we use a spring-based layout to highlight the modularity of the hub genes across modules
g1 <- graph.adjacency(as.matrix(adj),mode="undirected",weighted=TRUE,diag=FALSE)  
layoutFR <- layout.fruchterman.reingold(g1,dim=2) ## Spring based layout
l <- layout.fruchterman.reingold(g1) * scaleFactor
hs <- hub_score(g1, weights=NA)$vector

as <- authority_score(g1, weights=NA)$vector

## Make an adjacency matrix and set some arbitrary cutoffs and parameters for plotting the network

adjMat <- bicor(datExprYeast1[keepcol,])
maxsize <- min(25,nrow(adjMat))
metsize <- min(10,nrow(adjMat))
adjMat[adjMat<0.2] <- 0
kME <- apply(adj,2,sum)
cutoff <- sort(kME,decreasing=TRUE)[maxsize]
cutoff2 <- sort(kME,decreasing=TRUE)[metsize]
tophubs <- kME>=cutoff
metgenes <- names(kME)[kME>=cutoff2]
keepgenes <- c(keepgenes,metgenes)
adjMat <- adjMat[tophubs,tophubs]
numcors <- min(200,(maxsize^2-maxsize)/2)
topcors <- sort(as.numeric(adjMat),decreasing=TRUE)[numcors]
adjMat[adjMat<=topcors] <- 0

## Create a network object to store the data, plot these networks in a circle plot
library(igraph)
#converting network to igraph object
library(intergraph)
class(g1)
#g <- asIgraph(dynamicMods)

maxsize <- min(25,nrow(adj))
metsize <- min(10,nrow(adj))
adj[adj<0.2] <- 0
gA <- graph.adjacency(as.matrix(adj[1:5,1:5]),mode="undirected",weighted=TRUE,diag=FALSE) ## Top 5 in the center
gB <- graph.adjacency(as.matrix(adj[6:maxsize,6:maxsize]),mode="undirected",weighted=TRUE,diag=FALSE) ## Additioanl genes on the periphery
layoutCircle <- rbind(layout.circle(gA),layout.circle(gB)) ## Construct layout
windows()
plot(layoutCircle)

## Here we use a spring-based layout to highlight the modularity of the hub genes across modules
g1 <- graph.adjacency(as.matrix(adj),mode="undirected",weighted=TRUE,diag=FALSE)  
layoutFR <- layout.fruchterman.reingold(g1,dim=2) ## Spring based layout
## Plot the network for this module, specifying plotting parameters
#g <- graph.adjacency(A)
V(g1)$name <- rownames(adj)
V(g1)$name
windows()
plot.igraph(g1,layout=layoutFR, vertex.color="green", vertex.label=NA, 
            vertex.label.font=10, vertex.label.cex=.6)

par(mar=c(2,2,2,1))

###Visualize network for blue and brown modules
# Cytoscape
# select modules
modules = c("blue","brown")
probes = names(datExprYeast1)
modGenes=names(datExprYeast1)[modules = c("blue","brown")]
# Select module probes
dynamicColors = labels2colors(dynamicMods)
modulecolorsYCC=dynamicColors
inModule=is.finite(match(modulecolorsYCC,modules))
modProbes=probes[inModule]
march1 = match(modProbes, names(datExprYeast1))
modGenes=GeneAnnotation$gene_symbol[match1]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile=paste("CytoEdge",paste(modules,collapse="-"),
                                              ".txt",sep=""),
                               nodeFile=paste("CytoNode",paste(modules,collapse="-"),
                                              ".txt",sep=""),
                               weighted = TRUE, threshold = 0.02,nodeNames=modProbes,
                               altNodeNames = modGenes, nodeAttr = modulecolorsYCC[inModule])

#Hierarchical clustering dendrogram and module definition
windows()
par(mfrow=c(2,1), mar=c(2,2,2,1))
plot(geneTree, main="Yeast Gene Expression Network", labels=F, xlab="", sub="");
hclustplot1(geneTree,dynamicMods, title1="Colored by module membership")

collect_garbage()

cor.Y.rest <- bicor(datExprYeast1, maxPOutliers = 0.1)
adj2 <- adjacency.fromSimilarity(cor.Y.rest, type = "signed hybrid", power = softPower)
diss1=1-TOMsimilarity(adjMat =adj2, 
                      TOMType = "signed",
                      TOMDenom = "mean")

colnames(diss1) =rownames(diss1) =rownames(Y2)[restGenes]
hier1=hclust(as.dist(diss1), method="average" )
windows()
cex1=1
par(mfrow=c(1,2))
plotDendroAndColors(hier1, 
                    dynamicMods, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")


windows()
par(mfrow=c(1,2))
plotDendroAndColors(hier1, dynamicMods, "", autoColorHeight = FALSE,
                    colorHeight = 0.2,
                    ylab = "", xlab = "", sub = "", dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, cex.rowText = 1.3, cex.colorLabels = 1.2,
                    guideHang = 0.05, 
                    main = "Yeast GE module colors")

beta1=5
Connectivity=softConnectivity(datExpr,power=beta1)
ConnectivityCut = 1264 # number of most connected genes that will be considered
ConnectivityRank = rank(-Connectivity)
restConnectivity = ConnectivityRank <= ConnectivityCut
# thus our module detection uses the following number of genes
sum(restConnectivity)

ADJrest = adjacency(datExpr[,restConnectivity], power=beta1)
# The following code computes the topological overlap matrix based on the
# adjacency matrix.
dissTOM=TOMdist(ADJrest)
gc()

colorlevel1=levels(factor(dynamicMods))
summary(dynamicMods)
GeneSignificance=datSummary$essentiality[restConnectivity]
par(mfrow=c(2,1), mar=c(2,4,2,2))
windows()
plot(hierTOM, main="", labels=F, xlab="", sub="");
#plotColorUnderTree(hierTOM,colors=data.frame(dynamicMods,Essentiality = GeneSignificance))

####Heatmap for all modules
# for visualizing the network. Here we chose 2 scaling dimensions
windows()
cmd1=cmdscale(as.dist(diss1),2)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(dynamicMods),  main="MDS plot", xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

windows()
TOMplot1(diss1 , hier1, dynamicMods)

# Calculate eigengenes 
MEList = moduleEigengenes(datExprYeast1, colors = dynamicMods) 
MEs = MEList$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDiss = 1-cor(MEs); 
# Cluster module eigengenes 
METree = hclust(as.dist(MEDiss), method = "average"); 
# Plot the result 
windows() 
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

##Module membership
geneModuleMembership2 = signedKME(datExprYeast1, MEs)
colorsA2 = names(table(dynamicMods))
colnames(geneModuleMembership2)=paste("PC",colorsA2,".cor",sep="")
MMPvalue2 = corPvalueStudent(as.matrix(geneModuleMembership2),dim(datExprYeast1)[[2]]);
colnames(MMPvalue1)=paste("PC",colorsA2,".pval",sep="");
Gene = rownames(datExprYeast1)
kMEtable2 = cbind(Gene,Gene,dynamicMods)
for (i in 1:length(colorsA2))
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i],
                    MMPvalue2[,i])
colnames(kMEtable2)=c("PSID","Gene","Module",sort(c(colnames(
  geneModuleMembership2), colnames(MMPvalue2))))

##Comparing the 2 networks 
Gene = rownames(datExprADJ1[restDegreeADJ])
topGenesKME = NULL
for (c in 1:length(colorsA1)){
  kMErank1 = rank(-geneModuleMembership1[,c])
  kMErank2 = rank(-geneModuleMembership2[,c])
  maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = names(colorh1)
topGenesKME

#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge 
MEDissThres = 0.25 
# Call an automatic merging function 
merge = mergeCloseModules(datExprYeast1, dynamicMods, cutHeight = MEDissThres, verbose = 3) 
# The merged module colors 
mergedColors = merge$colors;
# Eigengenes of the new merged modules: 
mergedMEs = merge$newMEs;
#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
windows()
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6) 
plotDendroAndColors(hierTOM, cbind(dynamicMods, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
# Rename to moduleColors 
CCmoduleColors = mergedColors 
# Construct numerical labels corresponding to the colors 
colorOrder = c("grey", standardColors(50)); 
CCmoduleLabels = match(CCmoduleColors, colorOrder)-1;
table(CCmoduleLabels)
moduleColors2 = labels2colors(dynamicMods)
table(moduleColors2)
CCMEs = mergedMEs; 
geneTree = merge$dendrograms[[1]];

####Module differential Analysis
library(MODA)
ResultFolder = 'ConditionSpecificModules' 
# where middle files are stored
CuttingCriterion =  'Density' 
# could be Density or Modularity
indicator1 =  'X' 
# indicator for data profile 1
indicator2 =  'Y' 
# indicator for data profile 2
specificTheta = 0.7
#threshold to define condition specific modules
conservedTheta = 0.7 
#threshold to define conserved modules
##modules detection for network 1
intModules1 <- WeightedModulePartitionDensity(datExprYeast,ResultFolder,indicator1,CuttingCriterion)
##  ..done.
##modules detection for network 2
intModules2 <- WeightedModulePartitionDensity(ADJ,ResultFolder,indicator2,CuttingCriterion)
JaccardMatrix <- comparemodulestwonets(ResultFolder,intModules1,intModules2,
                                       paste('/DenseModuleGene_',indicator1,sep=''),
                                       paste('/DenseModuleGene_',indicator2,sep=''))
CompareAllNets(ResultFolder,intModules1,indicator1,intModules2,
               indicator2,specificTheta,conservedTheta)


###Differential analysis
library(WGCNA) 
datExprB1 <- ADJ1
datExprB2 <- dat0
probesI = NodeNames 
Yeastdata <- na.omit(dat0)
dim(Yeastdata)
datExprYeast1=as.data.frame(t(Yeastdata[,-c(1:7)]))
names(datExprYeast1)
probesA = names(datExprYeast1)
#Get gene symbols for Yeast PPI dataset
library(org.Sc.sgd.db)
Symbol <- mget(probesI, org.Sc.sgdGENENAME)
genesI = Symbol
###ID = unlist(mget(probesA,org.Sc.sgdGENENAME))
Symbol1 <- mget(probesA, org.Sc.sgdGENENAME)
genesA = Symbol1
datExprB1g = (collapseRows(datExprB1,genesI,probesI))[[1]] 
datExprB2g = (collapseRows(datExprB2,genesA,probesA))[[1]] 
#limit your analysis to genes/probes that are expressed in both data sets
commonProbesA = intersect (rownames(ADJ1),rownames(dat0))
datExprA1p = ADJ1[commonProbesA,] 
datExprA2p = dat0[commonProbesA,]
#commonGenesB = intersect (genesI,genesA) 
commonGenesB = intersect (rownames(datExprB1g),rownames(datExprB2g))
datExprB1g = datExprB1g[commonGenesB,] 
datExprB2g = datExprB2g[commonGenesB,] 
#Correlating general network properties
#Pick soft threshold
powers1=c(seq(1,10,by=1), seq(12,20, by=2))
RpowerTable=PickSoftThreshold(datExprYeast,powervector=powers1)[[2]]
collect_garbage()
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red")

softPower = 10 # (Read WGCNA tutorial to learn how to pick your power) 
rankExprA1= rank(rowMeans(datExprA1p)) 
rankExprA2= rank(rowMeans(datExprB2g)) 
randomcommongenes= sample(commonProbesA,677) 
rankConnA1= rank(softConnectivity(t(datExprA1p[randomcommongenes,]),type="signed",power=softPower)) 
rankConnA2= rank(softConnectivity(t(datExprA2p[randomcommongenes,]),type="signed",power=softPower)) 
rankExprB1= rank(rowMeans(datExprB1g)) 
rankExprB2= rank(rowMeans(datExprB2g)) 
randomcommongenes= sample(commonGenesB,677) 

# Now we find the whole network connectivity measures for each:
rankConnB1= rank(softConnectivity(t(datExprB1g[randomcommongenes,]),type="signed",power=softPower)) 
rankConnB2= rank(softConnectivity(t(datExprB2g[randomcommongenes,]),type="signed",power=softPower)) 
pdf("NetworkProperties.pdf", height=10, width=9) 
par(mfrow=c(2,2)) 
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (A1)",  
                   ylab="Ranked Expression (A2)") 
verboseScatterplot(rankConnA1,rankConnA2, xlab="Ranked Connectivity (A1)",  
                   ylab="Ranked Connectivity (A2)") 
verboseScatterplot(rankExprB1,rankExprB2, xlab="Ranked Expression (B1)",  
                   ylab="Ranked Expression (B2)") 
verboseScatterplot(rankConnB1,rankConnB2, xlab="Ranked Connectivity (B1)",  
                   ylab="Ranked Connectivity (B2)") 
dev.off() 

#choose the top 677 most expressed probes in data set A1
#and then keep only 1 probe per gene
keepGenesExpr = rank(-rowMeans(datExprA1p))<=677 
datExprA1g    = datExprA1p[keepGenesExpr,] 
datExprA2g    = datExprA2p[keepGenesExpr,] 
keepGenesDups = (collapseRows(datExprA1g,genesI,probesI))[[2]] 
datExprA1g    = datExprA1g[keepGenesDups[,2],] 
datExprA2g    = datExprA2g[keepGenesDups[,2],] 
rownames(datExprA1g)<-rownames(datExprA2g)<-keepGenesDups[,1] 

#calculate all of the necessary values to run WGCNA
library(WGCNA)

cor.Y <- bicor(datExprYeast1, use = "pairwise.complete.obs", maxPOutliers = 0.1)
adj = adjacency.fromSimilarity(cor.Y, type = "signed hybrid", power = softPower)
TOM = TOMsimilarity(adj, TOMDenom = "mean", TOMType = "signed")
colnames(TOM) <- rownames(TOM) <- rownames(Yeastdata)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM),method="average");

##biweight midcorrelation
cor.Z <- bicor(t(datExprA1g), use = "pairwise.complete.obs", maxPOutliers = 0.1)
adjacencyA1 = adjacency.fromSimilarity(cor.Z, type="signed hybrid", power=softPower);
diag(adjacencyA1)=0 
dissTOMA1   = 1-TOMsimilarity(adjacencyA1, TOMType="signed") 
library(flashClust)
geneTreeA1  = flashClust(as.dist(dissTOMA1), method="average") 
cor.V <- bicor(t(datExprA2g), use = "pairwise.complete.obs", maxPOutliers = 0.1)
adjacencyA2 = adjacency.fromSimilarity(cor.V, type="signed hybrid", power=softPower); 
diag(adjacencyA2)=0 
dissTOMA2   = 1-TOMsimilarity(adjacencyA2, TOMType="signed") 
geneTreeA2  = flashClust(as.dist(dissTOMA2), method="average") 

adjacencyA1 = adjacency(t(datExprA1g),power=softPower,type="signed hybrid"); 
diag(adjacencyA1)=0 
dissTOMA1   = 1-TOMsimilarity(adjacencyA1, TOMType="signed") 
geneTreeA1  = flashClust(as.dist(dissTOMA1), method="average") 
adjacencyA2 = adjacency(t(datExprA2g),power=softPower,type="signed hybrid"); 
diag(adjacencyA2)=0 
dissTOMA2   = 1-TOMsimilarity(adjacencyA2, TOMType="signed") 
library(flashClust)
geneTreeA2  = flashClust(as.dist(dissTOMA2), method="average") 
save.image("tutorial.RData")  #  (Section will take ~5-15 minutes to run) 

# display the networks visually
pdf("dendrogram.pdf",height=6,width=16) 
par(mfrow=c(1,2)) 
plot(geneTreeA1,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (A1)", 
     labels=FALSE,hang=0.04); 
plot(geneTreeA2,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (A2)", 
     labels=FALSE,hang=0.04);  
dev.off() 

##Heatmap of gene expressions of modules 
ClusterSamples=hclust(dist(datExprYeast1),method="average")
allcolors=unique(dynamicMods)
#for(i in 1:length(allcolors)){
which.module=allcolors
windows()
plotMat(t(scale(datExprYeast1[ClusterSamples$order,][,dynamicMods==which.module ]) ),nrgcols=30,rlabels=T, clabels=T,rcols=which.module,
        title=paste("heatmap of the module expressions" ) )
# Now we plot the expression values of the module eigengene
barplot(PC1$PCturquoise[ClusterSamples$order],col=which.module, main="Expression of the turquoise module eigengene")



#######Overlappig Modules PPI vs CC 
ng = (moduleLabels!="grey")&(CCmoduleLabels!="grey");
labels1 <- moduleLabels[ng];
labels2 <- CCmoduleLabels[ng];
overlap = overlapTable(labels1, labels2);
# Prepare axis labels for the table plot
xLabels = spaste("C.", sort(unique(labels2)), " (", spaste("ME", labels2colors(sort(unique(labels2)))), ")");
yLabels = spaste("P.", sort(unique(labels1)), " (", spaste("ME", labels2colors(sort(unique(labels1)))), ")");
# Content of the table
textMat = spaste(overlap$countTable, "|", signif(overlap$pTable, 2));
mat = overlap$pTable;
mat[mat<1e-60] = 1e-60;
lmat = -log10(mat);
#mat[mat>50] = 50;
dim(textMat) = dim(mat);
# Open a graphics window or a pdf file for plotting
x11()
par(cex = 0.8)
par(mar = c(11, 11, 5, 1));
labeledHeatmap(lmat,
               xLabels = spaste("ME", labels2colors(sort(unique(labels2)))),
               yLabels = spaste("ME", labels2colors(sort(unique(labels1)))),
               xSymbols = xLabels,
               ySymbols = yLabels,
               textMatrix = textMat, cex.text = 0.8,
               cex.main = 1.4,
               colors = greenWhiteRed(100)[50:100],
               setStdMargins = FALSE,
               main = "Overlap of PPI and Cell Cycle modules")
#dev.off();

#######cluster coefficient versus connectivity
CC= clusterCoef(ADJrest)
gc()
# Now we plot cluster coefficient versus connectivity
# for all genes
windows()
par(mfrow=c(1,1))
plot(Connectivity[restConnectivity],CC,col=as.character(dynamicMods),xlab="Connectivity",ylab="Cluster Coefficient" ) 


PPI= clusterCoef(AdjMat1rest)
gc()
# Now we plot cluster coefficient versus connectivity
# for all genes
windows()
par(mfrow=c(1,1))
plot(Degree[restDegreeADJ],PPI,col=as.character(colorh1),xlab="Connectivity",ylab="Cluster Coefficient" ) 

######Annotation
library(XML)
library(annotate)
library("org.Sc.sgd.db")
library(yeast2.db)
YPPIprobesI = NodeNames
gg <- getGO(Brn, "org.Sc.sgd.db")
sum(is.na(gg))

datOutput=data.frame(ProbeID=YPPIprobesI, gg)
write.table(datOutput,"YPPINewResults.csv",row.names=F, sep=",")
write.table(datOutput,"YPPINewResults.txt",row.names=F, sep=",")


# Read in the probe annotation 
annot <- as.matrix(read.table("(file location)/GO.ann.yeast.28.03.13.3_300.txt"), fill = TRUE)
# Match probes in the data set to the probe IDs in the annotation file 
dim(annot)
rownames(annot)[1:10]
colnames(annot)[1:2000]
datExprA=as.data.frame(t(annot[,-c(1:8)]))
str(annot)#Checking the output
probes2annot = match(probes, names(datExprA))
results <- cbind(probes, probes2annot)
head(results)

####Functional Similarity
library(GOSemSim)
YJR070C=c("GO:0016491","GO:0006519","GO:0000226") 
YAL015C=c("GO:0016829","GO:0016788", "GO.0016835") 
mgoSim(go1, go2, measure="Wang", ont = "CC")

YLR150W=c("GO:0048468","GO:0012501","GO:0016265") 
YHR064C=c("GO:0031326","GO:0006417", "GO.0010608") 
mgoSim(go1, go2, measure="Wang", ont = "CC")


#################GO Enrichment Analysis Fig Plot PPI and YGE
library(clusterProfiler)
library(org.Sc.sgd.db)
Symbol <- mget(data, org.Sc.sgdGENENAME)
gene.df <- bitr(Symbol, fromType = "GENENAME", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Sc.sgd.db)
head(gene.df)
ego2 <- enrichGO(gene = gene.df$ENTREZID, OrgDb = org.Sc.sgd.db, keytype = 'ENTREZID', ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
#head(summary(ego2))
head(summary(ego2), n=3) 
dotplot(ego2, showCategory=30)
enrichMap(ego2, vertex.label.cex=0.8, layout=igraph::layout.kamada.kawai)
x11()
par(cex = 0.8)
par(mar = c(11, 11, 5, 1));
plotGOgraph(ego2) 
X11()
par(cex = 0.8)
par(mar = c(11, 11, 5, 1));
barplot(ego2, drop=TRUE, showCategory=12)


Yeastdata <- na.omit(dat0)
dim(Yeastdata)
datExprYeast1=as.data.frame(t(Yeastdata[,-c(1:7)]))
probesA = names(datExprYeast1)
Symbol1 <- mget(probesA, org.Sc.sgdGENENAME)
#genesA = Symbol1
gene.df <- bitr(Symbol1, fromType = "GENENAME", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Sc.sgd.db)
head(gene.df)
ego2 <- enrichGO(gene = gene.df$ENTREZID, OrgDb = org.Sc.sgd.db, keytype = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
head(summary(ego2), n=3)
#Simplify is a method that reduces redundancy of enriched GO terms
ego2 <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
x11()
par(cex = 0.8)
par(mar = c(11, 11, 5, 1));
dotplot(ego2, showCategory=30)

X11()
par(cex = 0.8)
par(mar = c(11, 11, 5, 1));
barplot(ego2, drop=TRUE, showCategory=12)

##Protein Function Prediction
library(EGAD)
datExprADJ1 <- as.data.frame(t(ADJ1))
dim(datExprADJ1)
PPI <- names(datExprADJ1)
####Select annotations associated with our dataset from annotation data
x <- subset(annot, rownames(annot)%in%PPI)
write.csv(x,"PPI.csv",row.names=TRUE)

x <- subset(annot, rownames(annot)%in%PPI)
PP <- rownames(x)
YPP <- subset(ADJ1, rownames(ADJ1)%in%PP, colnames(ADJ1)%in%PP)
n <- nrow(YPP);
# Improve performance by Extend network
ext_gene_network <- extend_network(YPP, max=300)
annotations_sub1 <- extend_network(YPP)

#####scores 
annot_sub1 <- filter_network(x,flag = 1, min=100, max=1000)
ext_gene_network <- extend_network(YPP, max=300)

GO_groups_voted <- run_GBA(x, annotations_sub)
# neighbor voting
nv_results <- run_GBA(ext_gene_network, x)
head(nv_results)
mf_optimal <- calculate_multifunc(x)
windows()
hist(mf_optimal)
optimal_list<- as.numeric(mf_optimal[,4])
windows()
hist(optimal_list)
mf_results <- auc_multifunc(x, optimal_list)
gba_auc_nv <- neighbor_voting(x, ext_gene_network, nFold=3, output="AUROC")
gba_pr_nv <- neighbor_voting(x, ext_gene_network, nFold=3, output="PR")
head(gba_auc_nv)
head(gba_pr_nv)
multifunc_assessment <- calculate_multifunc(x)
auc_mf <- auc_multifunc(x, multifunc_assessment[,4])
windows()
hist <- plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)
filt <- !is.na(gba_auc_nv[,1])
aucA <- gba_auc_nv[filt,1]
aucB <- gba_auc_nv[filt,3]
windows()
hist <- plot_distribution(aucA, xlab="AUROCs")
windows()
avgs <- plot_density_compare(aucA, aucB, xlab="AUROC")
windows()
plot_value_compare(aucA, aucB)
scores <- predictions(annotations_sub1, ext_gene_network)
windows()
roc <- plot_roc(scores, annotations_sub1)
plot_roc_overlay(scores, annotations_sub1)
roc <- get_roc(scores, annotations_sub1)
auroc <- get_auc(roc[,1], roc[,2])
print(auroc)

##Evaluation
library(caret)
#table(scores)
xtab = table(factor(scores, levels=min(x):max(x)), factor(x, levels=min(x):max(x)))
# load Caret package for computing Confusion matrix
confusionMatrix(xtab)
(xtab)
# Precision: tp/(tp+fp):
xtab[1,1]/sum(xtab[1,1:2])
# Recall: tp/(tp + fn):
xtab[1,1]/sum(xtab[1:2,1])
# F-Score: 2 * precision * recall /(precision + recall):
2 * 1 * 1 / (1 + 1)

####Accuracy 
library(AUC)
accuracy(scores, x)
#### bar plots of Hub scores 
library(ggplot2)

df <- data.frame(Genes = c("ADK1","GAR1", "CDC5", "HHO1", "RAD27", "LSP1", "RPA49", "HSP104"), Hub_Score = c(0.956, 0.9213, 0.8021, 0.7826, 0.771, 0.6742, 0.6061, 0.5365))
df$Genes <- factor(df$Genes, levels = df$Genes[order(df$Hub_Score)])
df$Genes  # notice the changed order of factor levels
windows()
p <- ggplot(df, aes(x = Genes, y = Hub_Score, family="Times New Roman", width = 0.5)) +
  labs(x = "Genes") + 
  geom_bar(position = "dodge", stat='identity')  +
  theme(aspect.ratio = .1) + 
  scale_fill_grey() + theme_classic() +
  theme_bw(base_size = 12) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1.1,by=0.2))
p + coord_flip(ylim = c(0.5,1))

geom_col()


####Comparison of F-Score, Recall and Precision bar plots 
Unit <- c("Neighbor-Voting", "GBA", "CNPFP") 
Accuracy <- c(0.8879, 0.9386, 0.9710)
AUROC <- c(0.7241, 0.8502, 0.9862)
F_Score <- c(0.9406, 0.9683, 0.9691)
variable <- c(Accuracy, Precision, Recall, F_Score)
type <- c(rep("Accuracy", 3), rep("AUROC", 3), rep("F_Score", 3))
df <- data.frame(Unit, Accuracy, AUROC, F_Score)
head(df)
library(ggplot2)
require(tidyr)
library(reshape2)
df1 <- data.frame(Accuracy, AUROC, F_Score, unit)
df2 <- melt(df)
head(df2)
df.long <- gather(df, Unit)
colnames(df2) <- c("Unit", "variable", "value")
windows()
ggplot(data = df2, aes(x = Unit, y = value, fill = type, width = 0.50)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1.1,by=0.1))+
  theme(aspect.ratio = .1) + 
  scale_fill_grey() + theme_classic() +
  geom_col(position = position_dodge()) 
