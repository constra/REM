workwd <- "/Users/sun/Documents/Work_place/R/REM/"

### load data form Weissman 2003 expression data ###
Sc.exp.Weissman <- read.csv(paste(workwd,"DATA/Sc.exp.Weissman.csv",sep=""),na.strings="-", stringsAsFactor = F)
Sc.exp.Weissman[,"Protein.Molecules.Cell"] <- as.numeric(Sc.exp.Weissman[,"Protein.Molecules.Cell"])
str(Sc.exp.Weissman)

save(Sc.exp.Weissman,file = paste(workwd,"RData/Sc.exp.Weissman.RData",sep=""))
