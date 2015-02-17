### Fitness data from Qian et. al. 2012
fit.tab <- read.csv("DATA//Fitness_Qian2012.csv")
head(fit.tab)

## Find a gene
library(org.Sc.sgd.db)
ff <- function( ### ff stands for "fast filter"
    gn=NULL, ### gene name
    dat=NULL, ### name of the dataset
    col=NULL, ### which column are you looking for
    conv=TRUE, ### convert gene name
    target="ENSEMBL",src="GENENAME" ### convert gene name to ensembl
    ){
    if (conv){
        gs <<- select(org.Sc.sgd.db,keys=gn,columns=target,keytype=src)
        gs <- gs[[target]]
    } else{
        gs <-gn
    }
    
    f <- dat[which(dat[,col]==gs),]
    return(f)
}

gene <- "KSP1"
ff(gn=gene,dat=fit.tab,col="ORF")

### get SGD phenotype data from the web

#url <- url("http://downloads.yeastgenome.org/curation/literature/phenotype_data.tab")
url <- "http://downloads.yeastgenome.org/curation/literature/phenotype_data.tab"
phenodata <- readLines(url)
