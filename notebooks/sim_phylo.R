# Create M.bovis phylogeny from outbreak simulation
# D. Farrell May 2021

library(ape)
library(phangorn)
library('TransPhylo')
library('TreeDist')
setwd('~/gitprojects/pathogenie/notebooks/')

set.seed(2)
neg=100/365
off.r=3
w.shape=10
w.scale=0.1
pi=0.25
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2010,dateT=2014)
plot(simu)

ttree<-extractTTree(simu)
plot(ttree)
ptree<-extractPTree(simu)
plot(ptree)
p<-phyloFromPTree(ptree)
plot(p)
axisPhylo(backward = F)
write.tree(p,'sim.newick')

#load ref tree back in
reftree <- read.tree('sim.newick')
reftree <- root(reftree,'1')
plot(reftree)

#load snp and mlst trees
snptree <- read.tree('sim_results3/RAxML_bestTree.variants')
snptree <- drop.tip(snptree,'ref')
#snptree <- root(snptree,'1')
plot(snptree)
mlsttree <- read.tree('mlst3.newick')
#mlsttree <- root(mlsttree,'1')
plot(mlsttree)

#estimate mlst tree from dm
dm = read.table('dist_mlst.csv',sep=',',header=TRUE,row.names=1,check.names=FALSE)
dm
utree <- nj(as.matrix(dm))
plot(utree)

#compare trees
TreeDistance(reftree, mlsttree)
SharedPhylogeneticInfo(reftree, snptree)
VisualizeMatching(MutualClusteringInfo, reftree, mlsttree)
as.matrix(dm)

infer <- function(tree) {
  #infer ttree
  dateT=2014
  w.shape=10
  w.scale=0.1
  t <- tree
  t <- multi2di(t)
  t$edge.length <- pmax(t$edge.length,1/365) 
  ptree<-ptreeFromPhylo(t,dateLastSample=2013.5)
  plot(ptree)
  res<-inferTTree(ptree,mcmcIterations=1000,w.shape=w.shape,w.scale=w.scale,dateT=dateT)
  
  plot(res)
  med=medTTree(res)
  plot(med)
  myttree=extractTTree(med)
  plot(myttree,type='detailed',w.shape,w.scale)
}

infer(snptree)
