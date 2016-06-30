\name{RUN_GSAA}
\alias{RUN_GSAA}
\title{ RUN_GSAA }
\description{
  Master Function for running GSAA.}

\arguments{
    \item{input}{matrix with gene names as column names, group names as row names, and strengths as entries}
    \item{ensembl_dataset}{Which biomart will be used.  Species Specific.}
    \item{mart}{which specific dataset within the biomart will be used.  Platform Specific.}
    \item{GO_Member_Threshold}{Threshold for eliminating GO Terms with fewer than this many members from the incidence table (equivalent to a Minor Allele Frequency Filter in a GWAS)}
    \item{qval_threshold}{FDR threshold for significance.  Will take either 0-1 or 0-100 as inputs}
    }

\value{
    \item{Incidence_Table}{The Incidence_Table of genes and GO Terms (equivalent to a genotype file in GWAS)}
    \item{Sig_Results}{Significant GSAA results}
}

\usage{
  RUN_GSAA(input,ensembl_dataset="mmusculus_gene_ensembl",mart="illumina_mouseref_8_v2",GO_Member_Threshold=10,qval_threshold=.05)
}
\author{ Christoph Rau }
\note{
}
\seealso{ \code{\link{GSAA}},\code{\link{calculateIncidenceTable}}}