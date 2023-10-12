# code to convert vcf to distance matrix for Splitstree 
# Jenna McCullough & Devon DeRaad 

library(vcfR)
library(adegenet)
library(StAMPP)


# SACER 63 

# i'm giving the code for just the sacer 63 dataset. you can add in the subsets
# if you'd like to use them instead. 
v<-read.vcfR("sacer.63-filtered.thinned.vcf.gz")

#prepare input matrix

#convert vcfR to genlight
gen<-vcfR2genlight(v)
#check sample names (remember the splitstree gui won't let you visualize more 
#   than 10 characters in a sample name, 
# and all sample names must be unique)
gen@ind.names
#edit your sample names to fit this criteria. If you need to change names for
# some reason, you can do something like this
#gen@ind.names<-gsub("D_hypoleucum","hy", gen@ind.names)
#gen@ind.names<-gsub("D_nigrilore","ny", gen@ind.names)
#double check sample names
#gen@ind.names

#assign sample names as populations (population assignments are a requirement for the stampp functions, but can be arbitrary here because we want a pairwise divergence matrix among all samples, not samples assigned populations)
pop(gen) <- gen@ind.names

#make pairwise divergence matrix among all samples
sample.div <- stamppNeisD(gen, pop = FALSE)

#export for splitstree
stamppPhylip(distance.mat=sample.div, file="sacer63-distanceMatrix.txt")

# now open up splitstree gui and open the file and it will generate automatically 

