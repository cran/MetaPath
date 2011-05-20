\name{plotMAPE}
\alias{plotMAPE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot MAPE outcomes
}
\description{
This function will plot two figures. The first figure is the Venn diagram to show the overlapped enriched pathways identified by MAPE_G, MAPE_P and MAPE_I. 
The second figure is the heatmap of the q-values of enriched pathways. 
}
\usage{
plotMAPE(MAPE.obj, cutoff, MAPE.method = c("MAPE_I", "MAPE_P", "MAPE_G"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAPE.obj}{
The output of MAPE.
}
  \item{cutoff}{
The q-value cutoff.
}
  \item{MAPE.method}{
The MAPE method of interest. 
}
}

\value{
A Venn diagram and a heatmap of q-values of enriched pathways will be plotted.
When plot the heatmap, if the MAPE.method is MAPE_I, it will plot the q-values of enriched pathways for each individual study and q-values computed by three MAPE methods.
if the MAPE.method is MAPE_P, it will plot the q-values of enriched pathways for each individual study and q-values computed by the MAPE_P method.
if the MAPE.method is MAPE_G, it will plot the q-values of enriched pathways for each individual study and q-values computed by the MAPE_G method.
}

\examples{
\dontrun{
plot.MAPE(MAPE.obj, cutoff=0.05, MAPE.method = "MAPE_I") 
}
}

