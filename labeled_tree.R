setwd("~/Documents/NESCent Project/Main conus files/Latest versions") #set working directory

library(ape) #required for pruning tree

conus_56spp=read.csv("conus_56spp_forpub.csv", row.names = 1) #attach full Conus dataset
attach(conus_56spp)


############Produce pruned Conus tree############
conetree <- read.nexus("conetree.nex")

remove <- conetree$tip.label[-which(conetree$tip.label %in% rownames(conus_56spp))] #find mismatched tips

pru.conetree <- drop.tip(conetree, tip=remove) #drop names not in conus_56spp

######SWITCH CSV FILE SO THAT ROW NAMES ARE IN TREE NUMBER ORDER#########
##attach species names and log SL max values to pruned tree
log.max.ed.tree <- conus_56spp$log.max.ed
names(log.max.ed.tree) <- pru.conetree$tip.label
log.max.ed.tree <- log.max.ed.tree[match(pru.conetree$tip.label, names(log.max.ed.tree))]

# Draw tree with various characteristics plotted on it
plot(ladderize(pru.conetree), type = "phylogram", show.tip.label =T, label.offset = 3.0, edge.width = 2, cex = 1.1)
#To show distribution of non-planktonic hatchlings
tiplabels(pch=15, cex=1.5*crawl.juv, col = "gray45") 
# To show distribution of egg diameters
tiplabels(pch=1, cex=0.006*max.egg.diam, lwd=2.25, col = "black") #used deepskyblue4 before
# To show distribution of values of residuals
tiplabels(pch=1, cex=(resid+1)*3, lwd=2.25, col = "black")
# To show distribution of # eggs per capsule
tiplabels(pch=1, cex=0.05*sqrt(max.eggs.cap), lwd=2.25, col = "black")
# To show distribution of predation types
tiplabels(pch=19, cex=pisc.pred, col = "aquamarine4")
tiplabels(pch=19, cex=verm.pred, col = "orchid4")
tiplabels(pch=19, cex=mol.pred, col = "tomato3")
legend("bottomleft", legend=c("fish prey", "worm prey", "mollusc prey"), fill=c("aquamarine4", "orchid4", "tomato3"), cex=0.9, bty="o")