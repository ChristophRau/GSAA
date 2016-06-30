\name{calculateIncidenceTable}
\alias{calculateIncidenceTable}
\title{ calculateIncidenceTable }
\description{
  Calculates a Gene/GO Term incidence Matrix (equivalent to a SNP Genotype file for a GWAS)}

\arguments{
  \item{input}{matrix with gene names as column names, group names as row names, and strengths as entries}
  \item{ensembl_dataset}{Which biomart will be used.  Species Specific.}
  \item{mart}{which specific dataset within the biomart will be used.  Platform Specific.}
  \item{GO_Member_Threshold}{Threshold for eliminating GO Terms with fewer than this many members from the incidence table (equivalent to a Minor Allele Frequency Filter in a GWAS)}
  }

\value{
  \item{Incidence_Table}{The Incidence_Table of genes and GO Terms (equivalent to a genotype file in GWAS)}
  }

\usage{
  calculateIncidenceTable(input,ensembl_dataset="mmusculus_gene_ensembl",mart="illumina_mouseref_8_v2",GO_Member_Threshold=10)
}
\author{ Christoph Rau }
\note{
}
\seealso{ \code{\link{GSAA}},\code{\link{matchkeep}}}