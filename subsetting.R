###Subsetting###

getwd()
setwd("/Users/sun/Documents/Work/Work_place/R")
load("REM/RData/interactome.all.RData")
colnames(interactome.all)

### restrict only one term ###
interactome.Hela=interactome.all[,c("EnsemblID","Name","HeLa")]
interactome.Hela.RNAbinding=interactome.Hela[which(interactome.Hela[,"HeLa"] %in% "RNAbinding"),]

### restrict several terms ###
interactome.any=interactome.all
interactome.any.bind=interactome.any[which(interactome.any[,c("HeLa"] %in% "RNAbinding"),]

interactome.all.ensembl=as.vector(interactome.all[,"EnsemblID"])
#add description#
Hela.ensembl=as.vector(interactome.Hela.RNAbinding[,"EnsemblID"])

library(biomaRt)
ensembl=useMart("ensembl")
listDatasets(ensembl)

#get description for Hela#
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
Hela.ensembl.desp=getBM(attributes=c('ensembl_gene_id','description'),
                        filters = 'ensembl_gene_id',
                        values = Hela.ensembl, mart = ensembl)
#get description for all#
all.ensembl.desp=getBM(attributes=c('ensembl_gene_id','external_gene_id','description'),
                       filters='ensembl_gene_id',
                       values=interactome.all.ensembl,
                       mart=ensembl)
head(all.ensembl.desp)

## append description from biomaRt##
for(i in 1:length(interactome.Hela.RNAbinding[,"EnsemblID"]))
{
    replace.desp=Hela.ensembl.desp[which(Hela.ensembl.desp[,"ensembl_gene_id"] %in% interactome.Hela.RNAbinding[i,"EnsemblID"]),"description"]
    if(any(Hela.ensembl.desp[,"ensembl_gene_id"] %in% interactome.Hela.RNAbinding[i,"EnsemblID"]))
    {interactome.Hela.RNAbinding[i,"description"] = replace.desp
    }else{
        interactome.Hela.RNAbinding[i,"description"]=NA
    }
}

## write output table ##
write.table(interactome.Hela.RNAbinding,file="REM/OUTPUTS/interactome.Hela.RNAbinding.csv",sep="\t")
