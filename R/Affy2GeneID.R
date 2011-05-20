Affy2GeneID <-
function(madata,geneID='ENTREZID',annotation){

## load annotation data package
require('AnnotationDbi')

chip=strsplit(annotation,'.db')[[1]]
probe2gene =mget(featureNames(madata),get(paste(chip,geneID,sep='')))

# Remove probe identifiers that do not map to any ENTREZID
probe2gene <- probe2gene[!is.na(probe2gene)]
probe2gene <- probe2gene[!is.null(probe2gene)]


gene.uni=unique(unlist(probe2gene))
exprs.uni.IQR=matrix(NA,length(gene.uni),ncol(madata))
rownames(exprs.uni.IQR)=gene.uni

count=unlist(lapply(probe2gene,length))
rm.idx=which(count>1)
if(length(rm.idx)>1){
	probe2gene=probe2gene[-rm.idx]
}


probe2gene.mtx=matrix(NA,length(probe2gene),2)

probe2gene.mtx[,1]=names(unlist(probe2gene))
probe2gene.mtx[,2]=unlist(probe2gene)
madata.mtx=exprs(madata)[probe2gene.mtx[,1],]

multi.idx=rep(NA,length(gene.uni))
names(multi.idx)=gene.uni


for(t1 in 1:length(gene.uni)){

	expr.idx=madata.mtx[which(probe2gene.mtx[,2]==gene.uni[t1]),]
	if(is.null(nrow(expr.idx))){
	exprs.uni.IQR[t1,]=expr.idx
	multi.idx[t1]=0
	} else {
	id.IQR.idx=which.max(apply(expr.idx,1,IQR))
	exprs.uni.IQR[t1,]=expr.idx[id.IQR.idx,]

	
}
}


madata.IQR<- new("ExpressionSet", exprs = exprs.uni.IQR,phenoData=phenoData(madata) )

if(geneID==toupper('ENTREZID')) {
	featureNames(madata.IQR)=toupper(featureNames(madata.IQR))
}

out=madata.IQR 

return(out)
}
