#####GSAA#####
#
# This function calculates the associations, determines FDRs and returns the final results
# INPUTS:
#     input:  Same as above.
#     GoInfo: incidence matrix from above
#     qval_threshold:  FDR threshold.  Will take either 0-1 or 0-100 values and adjust accordingly
# OUTPUTS:  table of all significant GO Associations.  It will also generate a file ("GSAA_All_Pvalues.csv") with all pvalues.
#
###################

GSAA <- function(input, GoInfo,qval_threshold=.05){
  outfile<-file("GSAA_All_Pvalues.csv","w")
  cat("Module,GO Term,Pval",file=outfile,sep="\n")
  print(paste0("There Are ",nrow(input)," Groupings..."))
  for(i in 1:nrow(input)){
    print(paste0("Calculating Pvalues For Group ",i))
    bar=txtProgressBar(style=3)
    for(j in 1:nrow(GoInfo)){
      setTxtProgressBar(bar,j/nrow(GoInfo))
      pval=wilcox.test(as.numeric(input[i,])~as.numeric(GoInfo[j,]))$p.value #calculates the wilcoxan rank sum test score for the pvalue.
      outrow=c(rownames(input)[i],rownames(GoInfo)[j],pval)
      outrow=paste(outrow,collapse=",")
      cat(outrow,file=outfile,sep="\n")
    }
  }
  close(outfile)
  #reads it back in...
  pvals=read.csv("GSAA_All_Pvalues.csv")

  print("Calculating FDRs And Filtering")
  #Calculates the FDRs for the data.
  fdrs=fdrtool(pvals[,3],statistic="pvalue",plot=FALSE)
  pvals=cbind(pvals,fdrs$qval)

  if(qval_threshold>1){ qval_threshold=qval_threshold/100}

  p10=max(fdrs$pval[fdrs$qval<.1])
  p5=max(fdrs$pval[fdrs$qval<.05])
  p1=max(fdrs$pval[fdrs$qval<.01])
  pthresh=max(fdrs$pval[fdrs$qval<qval_threshold])
  print(paste0("Ten Percent FDR  P= ",p10))
  print(paste0("Five Percent FDR  P= ",p5))
  print(paste0("One Percent FDR  P= ",p1))
  print(paste0("Provided Percent FDR  P= ",pthresh))

  #Filters results to only keep those which surpass our FDR cutoff
  tokeep=fdrs$qval<qval_threshold
  print(paste0("Keeping ",sum(tokeep), " Elements With Pvalue Less Than ",pthresh))
  pvals=pvals[tokeep,]
  colnames(pvals)=c("Group","GO Term","P value","Q value (fdr)")

  print("Determining GO Term Names...")
  #Get the mappings of all GOIDs to GO Names (actual english...)
  xx <- as.list(GOTERM)

  #Map GOIDs to GO Names
  terms=c()
  for(i in 1:nrow(pvals)){
    temp=match(pvals[i,2],names(xx))
    terms=c(terms,Term(xx[[temp]]))
  }

  print("Sorting...")
  #Add to our original dataset, and finish outputting.
  outdata=cbind(pvals,terms)
  colnames(outdata)=c(colnames(pvals),"GO Term Name")
  outdata=outdata[ order(outdata[,1], outdata[,3]), ]
  print("Complete!")
  return(outdata)
}