testStrains <- as.factor(c("A224T", "N274Y", "N445K", "P421L", "S3C", "T41T","tRNA", "WT"))
strStrains <- as.list(levels(testStrains))
globe <- "Global"

for (testright in testStrains) {
  # assign( eval(testright), setNames(vector("list", length(testStrains)), testStrains)) #lists are treated as columns with separate rows
  holder <- setNames(vector("list", length(testStrains)), testStrains) # This place holder is dirty, but I don't know how to get around it.
  for (testleft in testStrains){
    if (testleft!=testright){
      deseq_results_strain <- paste(testleft,"_comparedTo_", testright) #"Running Deseq"
      holder[testleft] <- deseq_results_strain



      # holder[testright] <- deseq_results_strain
      # assign( base::get
    }
    assign( eval(testright), holder)
    #
    # assign(, holder) #Store in named row
  }
}
tableList <- lapply(strStrains, get) #Create a list of the 8 lists by strain name
names(tableList) <- strStrains
allByAll <- do.call(cbind, tableList) #working on cbind: it creates a matrix, not a data frame. Could be a problem.