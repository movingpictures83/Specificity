library(SINCERA)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
	## first iteration

sc <- readRDS(paste(pfix, parameters["genes", 2], sep="/"))
#In the first iteration, we selected genes for clustering using the following criteria
# genes expressed in at least 10 cells in at least 6 samples
min.samples <- as.integer(parameters["minsamples",2])
min.expression <- as.integer(parameters["minexpression",2])
min.cells <- as.integer(parameters["mincells",2])
specifity.threshold <- as.numeric(parameters["threshold",2])
obj <- prefilterGenes(sc, pergroup=TRUE, min.expression=min.expression, min.cells=min.cells, min.samples=min.samples)
# genes with at least 0.7 specificity in at least 6 samples
obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=specifity.threshold)

# set the selected genes for clustering
sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))

saveRDS(sc, outputfile)
}
