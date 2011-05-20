\name{MAPE}
\alias{MAPE}
\title{
   perform the Meta-Analysis for Pathway Enrichment (MAPE) analysis.
}
\description{
Description: This is the major function in the MetaPath package to implement the Meta-analysis for Pathway Enrichment (MAPE) methods introduced by Shen and Tseng (2010). The function automatically performs MAPE_G (integrating multiple studies at gene level), MAPE_P (integrating multiple studies at pathway level) and MAPE_I (a hybrid method integrating MAEP_G and MAPE_P methods).

In the simulation and real data analyses in the paper, MAPE_G and MAPE_P have complementary advantages and detection power depending on the data structure. In general, the integrative form of MAPE_I is recommended to use. In the case that MAPE_G (or MAPE_P) detects almost none pathway, the integrative MAPE_I does not improve performance and MAPE_P (or MAPE_G) should be used.
}
\usage{
MAPE(arraydata,pathway.DB,resp.type=c('twoclass','multiclass','continuous','survival'),stat=c('maxP','minP','rth','Fisher'),rth.value=NULL,
	permutation=c('sample','gene'),nperm=500,size.min=15,size.max=500,knn.neighbors=10,qvalue.cal=c('permute','estimate'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arraydata}{
  The arraydata is a list of microarray data sets. 
  Each microarray data set can be either an ExpressionSet or a list.
  If the microarray data set is a list, then it includes five elements as follows:
 	 1)x--exprs data
	 2)y-- the phenotype of interests
	 3)z-- censoring.status if applicable. 1 stands for the event occurred and 0 stands for censored.
	 4)geneid
	 5)samplename
  If the microarray data set is in an ExpressionSet format, the users need to 
   1) store the phenotype of interests in the slot 'label'.
   2) store the censor data is the slot 'censoring.status' if applicable

}
  \item{pathway.DB}{
 The pathway database in a GeneSetCollection format defined by GSEABase. 
 The pathway database can be downloaded from Broad institue (http://www.broadinstitute.org/gsea).
 PLEASE use the function 'getGmt' provided in the GSEABase package to load the pathway database.
}
  \item{resp.type}{
 The phenotype of interest.It is one of the four values: 'twoclass','multiclass','continuous', 'survival'.
}
  \item{stat}{
 The meta-analysis statistic to be used to combine two studies. It is one of the four values: 'minP','maxP','rth','Fisher'.
}
  \item{rth.value}{
 The value of the rth statistic if the meta-anlaysis statistic is 'rth'. For example, rth.value=0.6.
}
  \item{permutation}{
 The options for using sample permutation or gene permutation when performing enrichment analysis. it is one of the two values: 'gene' and 'sample'. The default option is sample permutation.
}
  \item{nperm}{
 Number of permutations to be performed.
}
  \item{size.min}{
 The minimum size of pathways to be considered. The default value is 15.
}
  \item{size.max}{
 The maximum size of pathways to be considered. The default value is 500.
 }
   \item{ knn.neighbors}{
	Number of neighbors to be used in the knn imputation method(default=10)
 }
   \item{qvalue.cal}{
The method to calculate the q-values. The default method is to calcuate the q-values based on the permutation method. If qvalue.cal='estimate', the q-values were estimated based on the Storey's method. 
 }
 
}

\value{
The qvalue and pvalue of each pathway. 
}
\references{
Shen, Kui, and George C Tseng. Meta-analysis for pathway enrichment analysis when combining multiple microarray studies. Bioinformatics (Oxford, England) 26, no. 10 (April 2010): 1316-1323. doi:10.1093/bioinformatics/btq148. http://www.ncbi.nlm.nih.gov/pubmed/20410053.
}
\author{
Kui Shen and George C Tseng.
}


\examples{
\dontrun{
library(MetaPath)
data(MAQC)
data(pathway.DB)
## Supposed we are interested in the ER related pathways, we first store the ER information in the slot 'label'. Then perform MAPE on this data set.
MAQC[[1]]$label=MAQC[[1]]$ER_status
MAQC[[2]]$label=MAQC[[2]]$ER_status
nperm=10 ## To save the computational time in this example, nperm was set to 10. The default value is 500.
MAPE.sample.obj=MAPE(arraydata=MAQC,pathway.DB=pathway.DB,resp.type="twoclass",stat='maxP',rth.value=NULL,nperm=nperm,permutation='gene',size.min=15,size.max=500)
cutoff=.1
subset(MAPE.sample.obj$qvalue,MAPE_I<=cutoff)
plotMAPE(MAPE.sample.obj,cutoff,MAPE.method='MAPE_I')
}
}
