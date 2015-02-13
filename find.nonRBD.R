library(yeast2.db)
library(org.Sc.sgd.db)

Sc.interactome.all <- mRNAinteractome$ENSEMBL
length(Sc.interactome.all)

Sc.exp.all <- Sc.exp.Weissman$ORF

no.interactome <- setdiff(Sc.exp.all,Sc.interactome.all)

Sc.exp.no.inter <- Sc.exp.Weissman[Sc.exp.Weissman$ORF %in% no.interactome,]

Sc.exp.no.inter.sort <- Sc.exp.no.inter[order(Sc.exp.no.inter$Protein.Molecules.Cell,decreasing = T),]

#example
head(keys(org.Sc.sgd.db,keytype="ENSEMBL"))
gene_names=Sc.exp.no.inter.sort$ORF
Sc.noRBD <- select(org.Sc.sgd.db,gene_names,columns=c("GENENAME","DESCRIPTION"),keytype="ENSEMBL")

for (i in 1:nrow(Sc.noRBD)){
    Sc.noRBD$EXP[i] <- Sc.exp.Weissman$Protein.Molecules.Cell[which(Sc.exp.Weissman$ORF == Sc.noRBD$ENSEMBL[i])]
}

write.table(Sc.noRBD,"OUTPUTS/Sc.noRBD.ens2exp.csv",sep=",")
