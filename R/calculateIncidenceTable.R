#####CalculateIncidenceTable#####
#
# This function creates an incidence table (similar to a genotype file) mapping genes to GO terms
# INPUTS:
#     input:  matrix with gene_ids as column names, groupings as row names, and weights as elements
#     ensembl_dataset:  Which biomart will be used.  Needs to be species specific!
#     mart:  which specific dataset within the biomart will be used.  Needs to be platform specific!
#     GO_Member_Threshold:  Will eliminate any GO term with fewer than this many members from the incidence table (equivalent of a minor allele frequency filter in GWAS)
#
# OUTPUTS:  nxm Incidence table, where n is the number of genes and m is the number of GO terms which pass the threshold filter.
#
###################
calculateIncidenceTable <- function(input, ensembl_dataset="mmusculus_gene_ensembl",mart="illumina_mouseref_8_v2",GO_Member_Threshold=10) {

  print("Creating Master GO TERM List...")
  ensembl = useMart("ENSEMBL_MART_ENSEMBL", host = "https://uswest.ensembl.org")
  ensembl = useDataset(ensembl_dataset,mart=ensembl)

  #identifies which Ids are present in your data, and matches them to the dataset to get a set of mappings from ID to GO Term
  IDs=colnames(input)
  goids = getBM(attributes=c(mart,'go_id'), filters=mart, values=IDs, mart=ensembl)

  print("Filtering GO TERM List...")
  #Identifies which GO Terms meet our numerical requirements (no enrichments on a single gene!)
  temp=names(which(table(goids[,2])>=GO_Member_Threshold))
  goids=matrix(unlist(apply(goids,1,matchkeep,temp,2)),ncol=2,byrow=TRUE)


  #get the names of all the remaining GO Terms
  goid_names=names(table(goids[,2]))
  goid_names=goid_names[-1]

  print("Constructing Incidence Matrix...")
  #And construct a large lookup matrix to easily map genes to GO Terms (our 'SNP' Matrix)
  outdata=array(0,c(length(goid_names),length(IDs)))
  bar=txtProgressBar(style=3)
  for(i in 1:length(goid_names)){
    setTxtProgressBar(bar,i/length(goid_names))
    #print(i/length(goid_names)*100)
    val=goid_names[i]
    temp=grep(paste0("^",val,"$"),goids[,2])
    ids=goids[temp,1]
    temp=match(ids,IDs)
    outdata[i,temp]=1
  }
  colnames(outdata)=IDs
  rownames(outdata)=goid_names
  print("Complete!")
  return(outdata)
}
