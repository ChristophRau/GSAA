#####MatchKeep#####
#
# This function looks at a row/column, determines if a given element has a match in a master list, and keeps the row/column if it does.
# INPUTS:  X: row/column vector.   names_to_keep: master list to match to.  element: which element to attempt to match with.
# OUTPUTS:  NULL if no match, x if match.
#
###################

matchkeep <- function(x,names_to_keep, element){
  if(!is.na(match(x[element],names_to_keep))){
    return(x)
  } else {
    return(NULL)
  }
}