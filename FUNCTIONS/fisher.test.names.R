#fisher test without names
fisher.test.without= function(testset,annolist,population,
                           sortby = "Pvalue", 										# "Pvalue" or "OddsRatio"
                           cutoff = 10^(-3),										# cutoff
                           min.size = 10,											# min size
                           max.size = 1000											# max size
)
{
  tab = matrix(NA,nrow=length(annolist),ncol=7)
  rownames(tab) = names(annolist)
  colnames(tab) = c("Pvalue","OddsRatio","ExpCount","Count","Size","Globalsize","Catsize")	
  for (i in 1:length(annolist)){
    annos = intersect(population,annolist[[i]])
    a = length(intersect(testset,annos))
    b = length(intersect(testset,setdiff(population,annos)))
    c = length(intersect(setdiff(population,testset),annos))
    d = length(intersect(setdiff(population,testset),setdiff(population,annos)))
    ftest = fisher.test(matrix(c(a,c,b,d),nrow = 2,dimnames = list(c("D","notD"),c("A","notA"))),alternative = "greater")
    tab[names(annolist)[i],] = c(ftest$p.value,ftest$estimate,length(testset)/length(population)*length(annos),length(intersect(testset,annos)),length(annos),
                                 (length(annos)*length(annos))/length(population),(length(intersect(testset,annos))*length(annos))/length(testset))
  }
  tab = tab[which(tab[,"Size"] >= min.size & tab[,"Size"] <= max.size),,drop = FALSE]
  tab = tab[which(tab[,sortby,drop = FALSE] <= cutoff),,drop = FALSE]
  tab = tab[order(tab[,sortby,drop = FALSE]),,drop = FALSE]
  return(tab)
}

### fisher.test.names ###
fisher.test.names = function(testset,annolist,population,
                                 sortby = "Pvalue",     								# "Pvalue" or "OddsRatio"
                                 cutoff = 10^(-3),										# cutoff
                                 min.size = 10,											# min size
                                 max.size = 1000											# max size
)
{
  tab = matrix(NA,nrow=length(annolist),ncol=8)
  rownames(tab) = names(annolist)
  colnames(tab) = c("Pvalue","OddsRatio","ExpCount","Count","Size","Globalsize","Catsize","GeneNames")
  for (i in 1:length(annolist)){
    annos = intersect(population,annolist[[i]])
    a = length(intersect(testset,annos))
    b = length(intersect(testset,setdiff(population,annos)))
    c = length(intersect(setdiff(population,testset),annos))
    d = length(intersect(setdiff(population,testset),setdiff(population,annos)))
    ftest = fisher.test(matrix(c(a,c,b,d),nrow = 2,dimnames = list(c("D","notD"),c("A","notA"))),alternative = "greater")
    tab[names(annolist)[i],] = c(ftest$p.value,ftest$estimate,length(testset)/length(population)*length(annos),length(intersect(testset,annos)),length(annos),
                                 (length(annos)*length(annos))/length(population),(length(intersect(testset,annos))*length(annos))/length(testset), gsub(", ",",",toString(intersect(testset,annos))))
  }
  
  tab = tab[which(as.numeric(tab[,"Size"]) >= min.size & as.numeric(tab[,"Size"]) <= max.size),,drop = FALSE]
  tab = tab[which(as.numeric(tab[,sortby]) <= cutoff),,drop = FALSE]
  tab = tab[order(as.numeric(tab[,sortby])),,drop = FALSE]
  
  #tab = tab[which(tab[,"Size"] >= min.size & tab[,"Size"] <= max.size),,drop = FALSE]
  #tab = tab[which(tab[,sortby,drop = FALSE] <= cutoff),,drop = FALSE]
  #tab = tab[order(tab[,sortby,drop = FALSE]),,drop = FALSE]
  
  return(tab)
}

####################
#fisher plot
####################
fisher.plot = function(tab,   							# DTA.fisher.test results
                           col = c("lightgrey","green","red","black","purple"), 	# vector of colors
                           alpha = NULL, 											# alpha value for color opacity
                           xlim = NULL, 											# xlimits, standard graphics parameter
                           ylim = NULL, 											# ylimits, standard graphics parameter
                           col.axis = "black",										# color of the axis
                           margin = 20, 											# left margin of the plot if rotate = TRUE
                           maxcat = 10, 											# maximal number of GO-terms to plot
                           rotate = TRUE		 									# should the plot be rotated 90 degrees clockwise
)
{
  if (dim(tab)[1] >= maxcat){tab = tab[1:maxcat,,drop = FALSE]} else {tab = tab[1:dim(tab)[1],,drop = FALSE]}
  terms = rownames(tab)
  tab = tab/tab[,"Size"]*100
  if (rotate){
    if (!is.null(xlim)){print("xlim argument will be ignored !")}
    if (is.null(ylim)){ylim = c(0,100)}
    col = convertcolor(col,alpha)
    if (!is.null(margin)){par(mar = c(5,margin,4,2)+0.1+1)}
    plot(1,col="white",ylim=c(0,2*length(terms)+1),xlim=ylim,ylab="",xaxt="n",xlab="",main="",axes=FALSE)
    abline(v=0,col="black")
    axis(2,at=seq(1.5,2*length(terms)-0.5,2),labels=terms,las=2,cex.axis=0.8,col="white",col.axis=col.axis)
    axis(1,at=seq(0,100,10),col.axis=col.axis)
    for (j in 1:length(terms)){rect(0,j+j-1,tab[j,"Size"],j+j,col=col[1])}
    for (j in 1:length(terms)){rect(0,j+j-1,tab[j,"Count"],j+j,col=col[2])}
    abline(v=tab[1,"ExpCount"],col=col[3],lwd=2)
    for (j in 1:length(terms)){lines(c(tab[j,"Globalsize"],tab[j,"Globalsize"]),c(j+j-1.2,j+j+0.2),lwd=2,lty=3,col=col[4])}
    for (j in 1:length(terms)){lines(c(tab[j,"Catsize"],tab[j,"Catsize"]),c(j+j-1.2,j+j+0.2),lwd=2,col=col[4])}
  } else {
    if (!is.null(xlim)){print("xlim argument will be ignored !")}
    if (is.null(ylim)){ylim = c(0,100)}
    col = convertcolor(col,alpha)
    plot(1,col="white",xlim=c(0,2*length(terms)+1),ylim=ylim,xlab="",xaxt="n",ylab="",main="",axes=FALSE)
    abline(h=0,col="black")
    axis(1,at=seq(1.5,2*length(terms)-0.5,2),labels=terms,las=2,cex.axis=0.8,col="white",hadj=0.6)
    axis(2,at=seq(0,100,10))
    for (j in 1:length(terms)){rect(j+j,0,j+j-1,tab[j,"Size"],col=col[1])}
    for (j in 1:length(terms)){rect(j+j,0,j+j-1,tab[j,"Count"],col=col[2])}
    abline(h=tab[1,"ExpCount"],col=col[3],lwd=2)
    for (j in 1:length(terms)){lines(c(j+j-1.2,j+j+0.2),c(tab[j,"Globalsize"],tab[j,"Globalsize"]),lwd=2,lty=3,col=col[4])}
    for (j in 1:length(terms)){lines(c(j+j-1.2,j+j+0.2),c(tab[j,"Catsize"],tab[j,"Catsize"]),lwd=2,col=col[4])}
  }
}



