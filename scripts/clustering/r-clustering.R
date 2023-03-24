setwd('/media/pgdrive/sharif/exosomians')

options(stringsAsFactors = F)
options(scipen = 999)

# install.packages('Rtsne')

library(data.table)
library(umap)
library(Rtsne)
library(ggplot2)
library(parallel)
library(coop) # Cosine similarity
library(dbscan)


evPredictionsFile = 'predictions/ExoCNN/ExoCNN.ev.extreme.90.unique.probabilities.csv'
evDna2VecEmbedsLatentLayerFile = 'dna2vec/ev.extreme.90.dna2vec.1to5.latent.csv'

#### Sequence, Edit distance, HClust ####
### Input: Sequneces
### Distance: Edit (Levenshtein) distance
### Clustering: HClust

evx = fread(evPredictionsFile)
evxEditd = adist(evx$seq)
rownames(evxEditd) = evx$id

evSeqsEditDistFile = 'data/rds/ExoCNN.ev.extreme.90.unique.seqs.edit-distance.rds'
# saveRDS(evxEditd, evSeqsEditDistFile)
# evxEditd = readRDS(evSeqsEditDistFile)


numOfClusters = 100

evxEditdHclust = hclust(as.dist(evxEditd))

evEditDistHClustPdfFile = sprintf('data/pdfs/ExoCNN.ev.extreme.90.unique.edit-distance.hclust.%d.pdf',
                                  numOfClusters)

hclust.plotPdf = function(hc, numOfClusters, size, file)
{
  pdf(width = size, height = size, file = file)
  plot(hc)
  rect.hclust(hc, k = numOfClusters)
  dev.off()
}

hclust.plotPdf(evxEditdHclust, numOfClusters, 1000, evEditDistHClustPdfFile)

evxEditdHclustClusters = data.table(id = evx$id,
                                    seq = evx$seq,
                                    cluster = cutree(evxEditdHclust, k = numOfClusters))

evEditDistHClustClustersFile = sprintf('data/csvs/ExoCNN.ev.extreme.90.unique.edit-distance.hclust.%d.clusters.csv',
                                       numOfClusters)
fwrite(evxEditdHclustClusters, evEditDistHClustClustersFile)


#### Dna2Vec, Cosine Distance, Dbscan ####
### Input: Dna2Vec Embeddings
### Distance: 1 - Cosine similarity 
### Clustering: Dbscan


evx = fread(evPredictionsFile)
evxD2v = fread(evDna2VecEmbedsLatentLayerFile)

numOfEmbedsPerSeq = 5
attr(evxD2v, 'id') = rep(evx$id, each = numOfEmbedsPerSeq)[!duplicated(evxD2v)]

evxD2v = unique(evxD2v)

#### Visulations via dimentionality reductions methods
evxD2vPca = prcomp(evxD2v)
evxD2vUmap = umap(as.matrix(evxD2v))
evxD2vTsne = Rtsne(evxD2v)

ggPca = data.table(x = evxD2vPca$x[, 1],
                   y = evxD2vPca$x[, 2])
ggplot(ggPca, aes(x, y)) + geom_point()

ggUmap = data.table(x = evxD2vUmap$layout[, 1],
                    y = evxD2vUmap$layout[, 2])
ggplot(ggUmap, aes(x, y)) + geom_point()

ggTsne = data.table(x = evxD2vTsne$Y[, 1],
                    y = evxD2vTsne$Y[, 2])
ggplot(ggTsne, aes(x, y)) + geom_point()



### Cosine similarity/distance
evxD2vCosd = cosine(t(as.matrix(evxD2v)))
evxD2vCosd = as.dist(1 - evxD2vCosd)


### Dbscan clustering
minValidNumOfClusterMembers = 30
eps = c(0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)

evxD2vCosdDbscan = mclapply(eps, function(eps)
  dbscan(evxD2vCosd, eps = eps, minPts = minValidNumOfClusterMembers),
  mc.cores = detectCores() - 2)

names(evxD2vCosdDbscan) = eps

evxD2vCosdDbscanClusters = lapply(evxD2vCosdDbscan, function(df) df$cluster)
lapply(evxD2vCosdDbscanClusters, table)

.PlotCluster = function(aCluster) {
  
  aClusterTable = as.data.table(table(aCluster))
  setorder(aClusterTable, N)
  
  aCluster = factor(aCluster, levels = aClusterTable$aCluster, ordered = T)
  ggUmap$cluster = aCluster
  ggTsne$cluster = aCluster
  
  setorder(ggUmap, -cluster)
  setorder(ggTsne, -cluster)
  
  ggplot(ggUmap, aes(x, y, color = cluster)) +
    geom_point() +
    scale_color_brewer(palette = 'Dark2') +
    ggtitle('UMAP', sprintf('Parameter = %s', aClusterName))
  
  ggplot(ggTsne, aes(x, y, color = cluster)) +
    geom_point() +
    scale_color_brewer(palette = 'Dark2') +
    ggtitle('TSNE', sprintf('Parameter = %s', aClusterName))
}


evD2vCosdDbscanPdfFile = 'data/pdfs/ExoCNN.ev.extreme.90.unique.cosine-dist.dbscan.pdf'
pdf(file = evD2vCosdDbscanPdfFile, width = plotSize, height = plotSize)

lapply(names(evxD2vCosdDbscanClusters), function(aClusterName) {
  aCluster = evxD2vCosdDbscanClusters[[aClusterName]]
  .PlotCluster(aCluster)
})

dev.off()

#### Dna2Vec, Cosine Distance, HClust ####
### Input: Dna2Vec Embeddings
### Distance: 1 - Cosine similarity 
### Clustering: Hclust

evxD2vCosdHclust = hclust(evxD2vCosd)

numOfClusters = c(5, 10, 30, 100, 500)
evxD2vCosdHclustClusters = mclapply(numOfClusters, function(k) cutree(evxD2vCosdHclust, k),
                                    mc.cores = detectCores() - 2)
lapply(evxD2vCosdHclustClusters, table)


evD2vCosdHclustPdfFile = 'data/pdfs/ExoCNN.ev.extreme.90.unique.cosine-dist.hclust.pdf'
pdf(file = evD2vCosdDbscanPdfFile, width = plotSize, height = plotSize)

lapply(names(evxD2vCosdDbscanClusters), function(aClusterName) {
  aCluster = evxD2vCosdDbscanClusters[[aClusterName]]
  .PlotCluster(aCluster)
})

dev.off()

#### Dna2Vec, Euclidean Distance, Dbscan ####
### Input: Dna2Vec Embeddings
### Distance: Euclidean
### Clustering: Dbscan

minValidNumOfClusterMembers = 30
eps = c(0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)

evxD2vEucdDbscan = mclapply(eps, function(eps)
  dbscan(evxD2v, eps = eps, minPts = minValidNumOfClusterMembers),
  mc.cores = detectCores() - 2)

names(evxD2vEucdDbscan) = eps

evxD2vEucdDbscanClusters = lapply(evxD2vEucdDbscan, function(df) df$cluster)
lapply(evxD2vEucdDbscanClusters, table)

evD2vEucdDbscanPdfFile = 'data/pdfs/ExoCNN.ev.extreme.90.unique.euclidean-dist.dbscan.pdf'
pdf(file = evD2vEucdDbscanPdfFile, width = plotSize, height = plotSize)

lapply(names(evxD2vEucdDbscanClusters), function(aClusterName) {
  aCluster = evxD2vEucdDbscanClusters[[aClusterName]]
  .PlotCluster(aCluster)
})

dev.off()


#### <NOT_IMPORTANT>
### Comparing cosine similarity and euclidean distance in
# two sets (containing 5 seqs which are generated from a similar seq) of sequences
# -> intra-set similarity is a little bit higher than inter-set one,
# and intra-set distance is a little bit lower than inter-set one,
# (at least most of the times)

i = 3
startIndex = 1 + (5 * i)
numOfReps = 5

evxD2vTest = evxD2v[1:(startIndex + (2*numOfReps) - 1), ]
evxD2vTest = t(as.matrix(evxD2vTest))

print('Cosine')
mean(as.numeric(cosine(evxD2vTest)[startIndex:(startIndex + (1*numOfReps) - 1),
                                   startIndex:(startIndex + (1*numOfReps) - 1)]))
mean(as.numeric(cosine(evxD2vTest)[(startIndex + (1*numOfReps)):(startIndex + (2*numOfReps) - 1),
                                   (startIndex + (1*numOfReps)):(startIndex + (2*numOfReps) - 1)]))
mean(as.numeric(cosine(evxD2vTest)[startIndex:(startIndex + (1*numOfReps) - 1),
                                   (startIndex + (1*numOfReps)):(startIndex + (2*numOfReps) - 1)]))

evxD2vTest = evxd2v[1:(startIndex + (2*numOfReps) - 1), ]

print('Euclidean')
mean(as.numeric(as.matrix(dist(evxD2vTest))[startIndex:(startIndex + (1*numOfReps) - 1),
                                            startIndex:(startIndex + (1*numOfReps) - 1)]))
mean(as.numeric(as.matrix(dist(evxD2vTest))[(startIndex + (1*numOfReps)):(startIndex + (2*numOfReps) - 1),
                                            (startIndex + (1*numOfReps)):(startIndex + (2*numOfReps) - 1)]))
mean(as.numeric(as.matrix(dist(evxD2vTest))[startIndex:(startIndex + (1*numOfReps) - 1),
                                            (startIndex + (1*numOfReps)):(startIndex + (2*numOfReps) - 1)]))
#### </>

