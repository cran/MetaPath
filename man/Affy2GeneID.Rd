\name{Affy2GeneID}
\alias{Affy2GeneID}

\title{
  Map Affymatrix probe IDs to other gene IDs.
}
\description{
  This function is used to map Affymatrix IDs to other gene IDs using Affymetrix annotation data (it requires the corresponding annotation data package). If a couple of probe IDs were mapped to same gene IDs, the probe ID which 
  has the largest IQR value over samples was selected to as the Representative of corresponded gene ID.
}
\usage{
Affy2GeneID(madata, geneID = "ENTREZID", annotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{madata}{
	Microarray data set in the ExpressionSet format.
}
  \item{geneID}{
	The geneIDs to be mapped. For example, 'ENTREZID' or 'SYMBOL'
}
  \item{annotation}{
 The annotation data pacakge provided by Bioconductor. For example, "hgu133a.db". 
 The annotation data packages are available at  http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData. The data packages can be installed by typing the following in an R command window:
  	source("http://bioconductor.org/biocLite.R")
	biocLite('hgu133a.db') ### if your chip is hgu133a.
The available map between different IDs can be found by typing 
    help(package=hgu133a.db)
}
}
\value{
The eSet with mapped gene IDs as featureNames.
}

\examples{
\dontrun{
require(hgu133a.db)
help(package=hgu133a.db)
eset=Affy2GeneID(madata=microarray.eset, geneID = "ENTREZID", annotation='hgu133a.db')
}
}





