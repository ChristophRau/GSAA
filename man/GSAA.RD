\name{GSAA}
\alias{GSAA}
\title{GSAA }
\description{
  Performs Gene Set Association Analysis, a weighted version of GSEA using GWAS ideas.}

\arguments{
  \item{input}{matrix with gene names as column names, group names as row names, and strengths as entries}
  \item{GOInfo}{Gene/GO Term Incidence Matrix (eqivalent to a SNP genotype file)}
  \item{qval_threshold}{FDR threshold for significance.  Will take either 0-1 or 0-100 as inputs}
}

\value{
  \item{Sig_Results} {Significant GSAA results}
}

\usage{
  GSAA(input,GOInfo,qval_threshold=.05)
}
\author{ Christoph Rau }
\note{
}
\seealso{ \code{\link{RUN_GSAA}},\code{\link{calculateIncidenceTable}}}