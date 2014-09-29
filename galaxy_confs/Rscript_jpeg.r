library(gplots)
library(lattice)

YlOrRd <- c("#001999", "#0022A3", "#0029B0", "#0032BB", "#0039C6", "#0050D2", "#0060E0", "#0070FF","#008AFF","#009AFF", "#00CAFF", "#FEE388","#FEE183","#FEDE7E","#FEDC7A","#FEDA78","#FED976", "#FECC69", "#FEC965", "#FEC561", "#FEC15B", "#FEBB57", "#FEB651", "#FEB24C", "#FEAB4B", "#FEAA4A", "#FEA849", "#FEA649", "#FEA448","#FEA247", "#FE9F46", "#FE9D45", "#FE9B45", "#FE9944", "#FE9843", "#FE9643", "#FE9442", "#FD9240", "#FD903F", "#FD8F3E", "#FD8E3D", "#FD8E3C", "#FD8D3C", "#FD8A3B", "#FD863A","#FD8339","#FD7F38","#FD7B37","#FD7936","#FD7635","#FD7234","#FC6E33","#FC6B32","#FC6931","#FC6530","#FC622F","#FC5E2E","#FE5A2D","#FC572C", "#FC4E2A", "#F94729", "#F63E27", "#F33725", "#F02E23", "#EC2721", "#EA2420", "#E9221F", "#E7201E", "#E61E1D", "#E31A1C", "#DE191D","#DC171D","#DA161E","#D9141E","#D7131F","#D51220","#D31120", "#CF0E21", "#CD0C21","#CB0B22","#CA0A22","#C80923", "#C60523","#C40024","#C20024","#C00025","#BD0025","#BB0026","#B80026","#B60026", "#B20026", "#AE0026", "#AB0026", "#A90026", "#A60026", "#A20026", "#9E0026", "#9B0026", "#990026", "#960026", "#920026", "#8E0026", "#8B0026", "#890026", "#860026", "#830026", "#800026")
YlOrRd.Lab <- colorRampPalette(YlOrRd, space = "Lab")
pmat <-read.delim("/home/gbotha/galaxy-dist/database/files/001/dataset_1429.dat",header=TRUE,row.names=1, skip=1)

x<-matrix(data=pmat[,1],ncol = 10, byrow = TRUE, dimnames=NULL)
use<-x
brks <- seq(-1,100)
labelY<-c("Vpr:1-10","Vpr:11-20","Vpr:21-30","Vpr:31-40","Vpr:41-50","Vpr:51-60","Vpr:61-70","Vpr:71-80","Vpr:81-90","Vpr:91-96")
jpeg(filename = "jpg_plot.jpg", quality=100)
heatmap.2(as.matrix(use), labRow=labelY, Rowv=FALSE, Colv=FALSE, dendrogram="none",breaks=brks, col=c("#000000", YlOrRd.Lab(length(brks)-2)), key=TRUE, symkey=FALSE, keysize=2,trace="none",density.info="none",na.color="#FFFFFF",margins=c(2,5),lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)), lwid=c(1,10,1),lhei=c(1.5,10,2.5),colsep=c(1,2,3,4,5,6,7,8,9),rowsep=c(1,2,3,4,5,6,7,8,9),sepcolor="#000000",sepwidth=c(0.001, 0.001),main="Heatmap of T0_Vpr",)
dev.off()
q()
n
