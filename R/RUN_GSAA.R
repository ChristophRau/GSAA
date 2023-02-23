#####RUN_GSAA#####
#
# Master function which calls the other functions in order to perform everything.
# INPUTS:
#     input:  Same as above.
#     ensembl_dataset:  Which biomart will be used.  Needs to be species specific!
#     mart:  which specific dataset within the biomart will be used.  Needs to be platform specific!
#     GO_Member_Threshold:  Will eliminate any GO term with fewer than this many members from the incidence table (equivalent of a minor allele frequency filter in GWAS)
#     qval_threshold:  FDR threshold.  Will take either 0-1 or 0-100 values and adjust accordingly
#
# OUTPUTS:  Returns Incidence Matrix AND Significant Results.
#
###################

RUN_GSAA <-function(input,ensembl_dataset="mmusculus_gene_ensembl",mart="illumina_mouseref_8",GO_Member_Threshold=10,qval_threshold=.05){
  library(fdrtool)
  library(GO.db)
  GOInfo=calculateIncidenceTable(input, ensembl_dataset=ensembl_dataset,mart=mart,GO_Member_Threshold=GO_Member_Threshold)
  SigResult=GSAA(input, GOInfo,qval_threshold=qval_threshold)
  newList = list("Incidence_Table"=GOInfo, "Sig_Results"=SigResult)
  return(newList)

}
