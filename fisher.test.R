######################
#Fisher test
#####################

prewd="/Users/sun/Documents/Work/Work_place/R/"
workdir=paste(prewd,"REM",sep="")
setwd(workdir)

annovector = list.files(paste(prewd,"ANNOTATION",sep=""))
for (i in annovector){load(paste(prewd,"ANNOTATION/",i,sep=""))}

source(paste("FUNCTIONS/","fisher.test.names.R",sep=""))

DTAfunctions = list.files(paste(prewd,"DTA",sep=""))
for (i in DTAfunctions){source(paste(prewd,"DTA/",i,sep=""))}

LSDfunctions = list.files(paste(prewd,"LSD",sep=""))
for (i in LSDfunctions){source(paste(prewd,"LSD/",i,sep=""))}

resolution = 2

table_of_interest="mRNAinteractome_evolution" #the file name of interest
table_import=read.csv2(paste("DATA/",table_of_interest,sep=""),head=T,sep="\t") #import the table of interest

### Analysing the file "mRNAinteractome_evolution"
testset=table_import[which(table_import$conserved..core..RBP == "core"),"ENSEMBL"]
testset=as.vector(testset) #gene list given by ensembl gene id
annolist=Sc.go.names.ensg
population=c(Sc.verified.ensg,Sc.uncharacterized.ensg)

list.fisher=fisher.test.names(testset,annolist,population)
list.fisher.nonames=fisher.test.without(testset,annolist,population)
write.table(list.fisher[,c("Pvalue","OddsRatio","ExpCount","Count","Size","GeneNames"),drop = FALSE],
            file = paste(workdir,"/OUTPUTS/",table_of_interest,"_go_","output",sep=""),na = "NA",append = FALSE, sep = "\t",col.names=NA,row.names=TRUE,quote = FALSE)

plotsfkt = function(){
    parfkt("default",1)
    if (dim(list.fisher)[1] > 0){
      fisher.plot(list.fisher.nonames,margin=20)
      xlab = expression(paste("Percentage"))
      ylab = expression(paste(""))
      main = paste("GO enrichment of",table_of_interest)
      #sub = paste("( ",length(repressed[[rownames(comparisons)[i]]])," of ",length(reliable)," genes tested )")
      mtextfkt("default",1,main,xlab,ylab)
    } else {emptyplot()}
}
plotit(filename = paste(workdir,"/PLOTS/",table_of_interest,"_plot",sep=""),sw = resolution*3,sh = resolution,sres = resolution,plotsfkt = plotsfkt,ww = 21,wh = 7,fileformat="pdf",saveit = TRUE,notinR = TRUE)

