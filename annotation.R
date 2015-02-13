######################
#Annotation
#####################

prewd="/Users/sun/Documents/Work/Work_place/R/"
workdir=paste(prewd,"REM",sep="")
setwd(workdir)

annovector = list.files(paste(prewd,"ANNOTATION",sep=""))
for (i in annovector){load(paste(prewd,"ANNOTATION/",i,sep=""))}

#source(paste("FUNCTIONS/","fisher.test.names.R",sep=""))

DTAfunctions = list.files(paste(prewd,"DTA",sep=""))
for (i in DTAfunctions){source(paste(prewd,"DTA/",i,sep=""))}

LSDfunctions = list.files(paste(prewd,"LSD",sep=""))
for (i in LSDfunctions){source(paste(prewd,"LSD/",i,sep=""))}
library(yeast2.db)
library(GO.db)

resolution = 2

table_of_interest="mRNAinteractome_evolution" #the file name of interest
table_import=read.csv2(paste("DATA/",table_of_interest,sep=""),head=T,sep="\t") #import the table of interest

### Analysing the file "mRNAinteractome_evolution"
testset=table_import[which(table_import$conserved..core..RBP == "core"),"ENSEMBL"]
testset=as.vector(testset) #gene list given by ensembl gene id

testset_anno=select(yeast2.db,testset,"GO","ENSEMBL")
select(GO.db,keys=testset_anno$GO,column="TERM",keytype="GOID")
