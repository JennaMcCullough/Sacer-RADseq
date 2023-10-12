# Code to make PCAs 
# Jenna McCullough 

# load packages
library(adegenet)
library(ade4)
library(vcfR)
library(scales)
library(parallel)
library(viridis)
library(wesanderson)
library(devtools)
library(ggplot2)

## Generating PCAs for all sacer samples and the subsets. 

## Sacer 63
# load in vcf, all of which are availiable on dryad, convert to genlight
sac63_vcf <- read.vcfR("01_sacer.63-filtered.thinned.vcf.gz")
### read in popmap that I used from STACKS, note that i turned it from a tab separated headerless file 
# to a comma separated file with headers "id" and "popmpap"
popmap<-read.csv("CARC-output/sacer63.csv")

sac63_gl <- vcfR2genlight(sac63_vcf)
pop(sac63_gl) <- c("vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "eximius", "eximius", "eximius", "eximius", "pealei", "juliae", "juliae", "juliae", "santoensis", "santoensis", "santoensis", "juliae", "santoensis", "ornatus", "solomonis", "solomonis", "solomonis", "solomonis", "solomonis", "amoenus", "amoenus", "amoenus", "amoenus", "amoenus", "pelewensis", "pelewensis", "farquhari", "farquhari", "farquhari", "farquhari", "farquhari", "leucopygius", "leucopygius", "leucopygius", "cimmamominus", "cimmamominus", "cimmamominus", "ruficollaris", "mauke", "mauke", "atiu", "gertrudae")
# runs PCA, retaining "nf" components
pca_sac63 <- glPca(sac63_gl, nf=20)

s.class(pca_sac63$scores[,1:2], pop(sac63_gl), col=c("brown1", "brown4", "chocolate4", "coral2", "aquamarine2", "aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "red", "orange", "yellow", "green", "dark green", "blue", "dark blue", "purple", "lavender", "palegreen3", "orchid4", "orchid2", "maroon3", "maroon4"), clab=1, cell=2.5)
pcaAllBarplot <- barplot(pca_sac63$eig/sum(pca_sac67$eig), main="eigenvalues")
head(pca_sac63$eig/sum(pca_sac63$eig))

# how i structure my pdf code is that i keep it at the bottom of each block of code so it generates
# a figure. I look at it, then i will run the PDF command,  plot it again, then run dev.off() 
# to make a PDF. i edit the color and style in illustrator to make them publication-ready. 
#pdf(file="sacer63.thinned.pca.pdf", width = 6, height = 6)
#dev.off()



####################
##### Sacer 45 #####
####################

# loads VCF using vcfR
sacer45_vcf <- read.vcfR("02_sacer45-filtered.thinned.vcf.gz")
# converts vcf to genlight
sacer45_vcf_gl <- vcfR2genlight(sacer45_vcf)
# assigns populations, in VCF order
pop(sacer45_vcf_gl) <- c("vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "eximius", "eximius", "eximius", "eximius", "pealei", "juliae", "juliae", "juliae", "santoensis", "santoensis", "santoensis", "santoensis", "santoensis", "ornatus", "solomonis", "solomonis", "solomonis", "solomonis", "solomonis", "amoenus", "amoenus", "amoenus", "amoenus", "amoenus")

# runs PCA, retaining "nf" components
pca_sacer45 <- glPca(sacer45_vcf_gl, nf=20)
# plots PCA with colors
s.class(pca_sacer45$scores[,1:2], pop(sacer45_vcf_gl), col=c("red", "orange", "yellow", "green", "blue", "purple", "black", "brown", "gray"), clab=1, cell=2.5)
# output barplot of eigenvalues to evaluate PCs
sacer45Barplot <- barplot(pca_sacer45$eig/sum(pca_sacer45$eig), main="eigenvalues")
head(pca_sacer45$eig/sum(pca_sacer45$eig))


pdf(file="sacer.45.thinned.pca.pdf", width = 6, height = 6)
dev.off()


####################
## PolySacer 26 ####
####################

# Sacer from polynesia with pealei 
# loads VCF using vcfR
fiji26_vcf <- read.vcfR("03_sacer26-filtered.thinned.vcf.gz")
# converts vcf to genlight
fiji26_vcf_gl <- vcfR2genlight(fiji26_vcf)
# assigns populations, in VCF order
pop(fiji26_vcf_gl) <- c("viti_VanuaLevu", "viti_VanuaLevu", "viti_Koro", "viti_Ovalau", "viti_Koro", "viti_Koro", "viti_Koro", "viti_Koro", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Kabara", "marin_Kabara", "marin_Kabara", "marin_OgeaDriki", "marin_OgeaLevu", "marin_OgeaLevu", "marin_OgeaLevu", "marin_Vanuvatu", "eximius", "eximius", "eximius", "eximius", "pealei")

# runs PCA, retaining "nf" components
pca_fiji26 <- glPca(fiji26_vcf_gl, nf=20)
# plots PCA with colors
s.class(pca_fiji26$scores[,1:2], pop(fiji26_vcf_gl), col=c("brown1","brown4","chocolate4","coral2","aquamarine2","aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "yellow"), clab=1, cell=2.5)
head(pca_fiji26$eig/sum(pca_fiji26$eig))
pdf(file="sacer-26-Fiji-peal.thinned.pca.pdf", width = 6, height = 6)
dev.off()
#

####################
## PolySacer 25 ####
####################

# Sacer from polynesia with no pealei 
# loads VCF using vcfR
fiji25_vcf <- read.vcfR("04_sacer25-filtered.thinned.vcf.gz")
# converts vcf to genlight
fiji25_vcf_gl <- vcfR2genlight(fiji25_vcf)
# assigns populations, in VCF order. You can change the population assignments easily by choosing to run any of these lines:
# this line will have them all with their own islands
pop(fiji25_vcf_gl) <- c("viti_VanuaLevu", "viti_VanuaLevu", "viti_Koro", "viti_Ovalau", "viti_Koro", "viti_Koro", "viti_Koro", "viti_Koro", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Kabara", "marin_Kabara", "marin_Kabara", "marin_OgeaDriki", "marin_OgeaLevu", "marin_OgeaLevu", "marin_OgeaLevu", "marin_Vanuvatu", "kadavu", "kadavu", "kadavu", "kadavu")
# this line will have them with their subspecies 
pop(fiji25_vcf_gl) <- c("vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "vitiensis", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "marinus", "eximius", "eximius", "eximius", "eximius")

# runs PCA, retaining "nf" components
pca_fiji25 <- glPca(fiji25_vcf_gl, nf=20)
# plots PCA with colors
s.class(pca_fiji25$scores[,1:2], pop(fiji25_vcf_gl), col=c("brown1","brown4","chocolate4","coral2","aquamarine2","aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "yellow"), clab=1, cell=2.5)
head(pca_fiji25$eig/sum(pca_fiji25$eig))
pdf(file="sacer-25-Fiji-peal.thinned.pca.pdf", width = 6, height = 6)
dev.off()
#

####################
##   Marinus 13 ####
####################
# loads VCF using vcfR
marinus_vcf <- read.vcfR("05_marinus13-filtered.thinned.vcf.gz")
# converts vcf to genlight
marinus_vcf_gl <- vcfR2genlight(marinus_vcf)
# assigns populations, in VCF order
pop(marinus_vcf_gl) <- c("marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Fulaga", "marin_Kabara", "marin_Kabara", "marin_Kabara", "marin_OgeaDriki", "marin_OgeaLevu", "marin_OgeaLevu", "marin_OgeaLevu", "marin_OgeaLevu", "marin_Vanuvatu")

# runs PCA, retaining "nf" components
pca_marinus <- glPca(marinus_vcf_gl, nf=20)
# plots PCA with colors
s.class(pca_marinus$scores[,1:2], pop(marinus_vcf_gl), col=c("brown1","brown4","chocolate4","coral2","aquamarine2","aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "yellow"), clab=1, cell=2.5)
head(pca_marinus$eig/sum(pca_marinus$eig))

#pdf(file="sacer-marinus_PC.thinned.pca.pdf", width = 6, height = 6)
#dev.off()
#

####################
## MelanSacer13 ####
####################
# just juliae, santoensis, and solomonis 

# load in vcf, convert to genlight
santo_vcf <- read.vcfR("07_santo13-filtered.thinned.vcf.gz")
# the difference between the different PCAs with and without the espiritu santo sample
# labelled with juliae or santoensis is just these lines. run the one you want to do: 
# below will have that sample with juliae 
pop(santo_gl) <- c("juliae", "juliae", "juliae", "santoensis", "santoensis", "santoensis", "juliae", "santoensis", "solomonis", "solomonis", "solomonis", "solomonis", "solomonis")
# below will have that sample with santoensis 
pop(santo_gl) <- c("juliae", "juliae", "juliae", "santoensis", "santoensis", "santoensis", "santoensis", "santoensis", "solomonis", "solomonis", "solomonis", "solomonis", "solomonis")
pca_santo <- glPca(santo_gl, nf=20)

s.class(pca_santo$scores[,1:2], pop(santo_gl), col=c("brown1", "brown4", "chocolate4", "coral2", "aquamarine2", "aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "red", "orange", "yellow", "green", "dark green", "blue", "dark blue", "purple", "lavender", "palegreen3", "orchid4", "orchid2", "maroon3", "maroon4"), clab=1, cell=2.5)
head(pca_santo$eig/sum(pca_santo$eig))

#pdf(file="MelanSacer13-as-juliae.pca.pdf", width = 6, height = 6)
dev.off()

####################
## MelanSacer13 ####
####################
#  juliae, santoensis,  solomonis and ornatus 

# load in vcf, convert to genlight
santo_vcf <- read.vcfR("06_santo.14-filtered.thinned.vcf.gz")
santo_gl <- vcfR2genlight(santo_vcf)
#below has the espiritu santo sample with juliae. 
pop(santo_gl) <- c("juliae", "juliae", "juliae", "santoensis", "santoensis", "santoensis", "juliae", "santoensis", "ornatus", "solomonis", "solomonis", "solomonis", "solomonis", "solomonis")
pca_santo <- glPca(santo_gl, nf=20)

s.class(pca_santo$scores[,1:2], pop(santo_gl), col=c("brown1", "brown4", "chocolate4", "coral2", "aquamarine2", "aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "red", "orange", "yellow", "green", "dark green", "blue", "dark blue", "purple", "lavender", "palegreen3", "orchid4", "orchid2", "maroon3", "maroon4"), clab=1, cell=2.5)
head(pca_santo$eig/sum(pca_santo$eig))

#pdf(file="MelanSacer14-as-juliae.pca.pdf", width = 6, height = 6)
#dev.off()

####################
## Micro10 ####
####################

# loads VCF using vcfR
micro_vcf <- read.vcfR("08_micro.10-filtered.thinned.vcf.gz")
# converts vcf to genlight
micro_vcf_gl <- vcfR2genlight(micro_vcf)
# assigns populations, in VCF order
pop(micro_vcf_gl) <- c("pelewensis","pelewensis","cimmamominus","cimmamominus","cimmamominus","ruficollaris","mauke","mauke","atiu","gert")

# runs PCA, retaining "nf" components
pca_micro <- glPca(micro_vcf_gl, nf=20)
# plots PCA with colors
s.class(pca_micro$scores[,1:2], pop(micro_vcf_gl),col=c("brown1","brown4","chocolate4","coral2","aquamarine2","aquamarine3", "cadetblue2", "cyan4", "cyan3", "darkorchid2", "yellow"), clab=1, cell=2.5)
head(pca_micro$eig/sum(pca_micro$eig))
#pdf(file="micro.thinned.pca.pdf", width = 6, height = 6)
#dev.off()
#
