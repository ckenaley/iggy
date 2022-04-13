library(geomorph)
library(tidyverse)
library(ggbiplot)
library(ggfortify)
library(fishtree)
library(ape)
library(treeio)

#https://academic.oup.com/iob/article/1/1/obz016/5526881#223239200
### https://gist.github.com/benmarwick/5799505

#devtools::install_github("ckenaley/geomorphcompanion")
library(geomorphcompanion)


meta <- read.csv("Dragonfish_inventory.csv")%>%
  mutate(cat=gsub(" ","",ID))

nts.list <- list.files(pattern = ".nts")
cat <- gsub("(MCZ\\d*).+","\\1",nts.list)

meta2 <- meta%>%
  select(cat,Genus,family,Diet)%>%unique()

meta2 <- tibble(cat=cat)%>%
  left_join(meta2)



lm.data <- array(dim=c(10,3,43))
for(i in 1:43){
  lm.data[,,i] <- readland.nts(nts.list[[i]])
  
}
dimnames(lm.data)[[3]] <- paste0(meta2$Genus)
lm.procrust <- gpagen(lm.data,ProcD = T)

 n.spec <- dim(lm.procrust$coords)[3]
 
 lm.list <- list()
for(i in 1:n.spec){
 lm.list[[i]] <-lm.procrust$coords[,,i]%>%as.tibble()%>%mutate(name=nts.list[[i]],pt=1:10)
}

lm.tbl <- do.call(rbind,lm.list)
lm.mean1 <- lm.tbl%>%
  mutate(cat=gsub("(MCZ\\d*).+","\\1",name))%>%
  left_join(meta2)%>%
  group_by(Genus,pt)%>%
  dplyr::summarize(x=mean(X),y=mean(Y),z=mean(Z))



lm.mean2 <- simplify2array(by(lm.mean1, lm.mean1$Genus, FUN = function(x) as.matrix(x[,-c(1:2)])))

head(lm.mean1)
lm.mean2[,,1]



gen <- dimnames(lm.mean2)[[3]]

diet <- tibble(Genus=gen)%>%
  right_join(meta2%>%select(Genus,Diet)%>%unique,copy = F)%>%select(Diet)%>%pull

stom.pca <- gm.prcomp(lm.mean2)
ggbiplot2(stom.pca,choices = c(1,2),labels = as.factor(gen))


##stom tree

kentree <- get.tree(treeio::read.beast("StomNucRYRhoBEAST50MilRun1-4meanhts.tree"))
kentree <- ladderize(kentree)
plot(kentree)


dups <- kentree$tip.label[duplicated(gsub("(.+)_.+","\\1",kentree$tip.label))]
kengen <- drop.tip(kentree,dups)

kengen$tip.label <- gsub("(.+)_.+","\\1",kengen$tip.label)
kengen <- drop.tip(kengen,"Thaleichthys")

drops <- setdiff(kengen$tip.label,meta2$Genus)
drops2 <- setdiff(gen,kengen$tip.label)
kengen <- drop.tip(kengen,c(drops,drops2))
plot(kengen)

inc <- which(gen%in%kengen$tip.label)
stom.Ppca <- gm.prcomp(lm.mean2[,,inc],phy=kengen,transform = F,GLS=T,align.to.phy = T)
plot(stom.Ppca)

gen2 <- dimnames(lm.mean2[,,inc])[[3]]

diet2 <- tibble(Genus=gen2)%>%
  right_join(meta2%>%select(Genus,Diet)%>%filter(Genus%in%gen2)%>%unique)%>%select(Diet)%>%pull%>%as.factor

fam2 <- tibble(Genus=gen2)%>%
  right_join(meta2%>%select(Genus,family)%>%filter(Genus%in%gen2)%>%unique)%>%select(family)%>%pull%>%as.factor

#phylo pca 1-2
pdf("PCA12.pdf")
plot(stom.Ppca,axis1 = 1,axis2=2,time.plot = T,tip.tex.cex=0.1,phylo.par = list(tip.labels = TRUE, tip.txt.cex = .75, edge.width = 2),bg=diet2,pch=22)
legend("topleft", pch=22, pt.bg = unique(diet2), legend = rev(levels(diet2)))
dev.off()


#phylo pca 1-3
pdf("PCA13.pdf")
plot(stom.Ppca,axis1 = 1,axis2=3,time.plot = T,tip.tex.cex=0.1,phylo.par = list(tip.labels = TRUE, tip.txt.cex = .75, edge.width = 2),bg=diet2,pch=22)
legend("topleft", pch=22, pt.bg = unique(diet2), legend = rev(levels(diet2)))

dev.off()

#morphological disprarity

gdf<- geomorph.data.frame(coords=lm.mean2[,,inc],diet=diet2,fam=fam2)

morphol.disparity(coords ~ 1, groups = NULL, data = gdf, 
                  iter = 999, print.progress = FALSE)

morphol.disparity(coords ~ 1, groups = ~diet, data = gdf, 
                  iter = 999, print.progress = FALSE)


morphol.disparity(coords ~ 1, groups = ~fam, data = gdf, 
                  iter = 999, print.progress = FALSE)
