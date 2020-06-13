library("DESeq2")
library("here")

read_table_counts <- function(filename){
	table_counts <- read.csv(filename, sep = "\t", row.names = 'locusTag')
	# Need to convert the table counts to integers because they may be float with zeros after the decimal point.	
	#table_counts[] <- lapply(table_counts, function(x) as.numeric(x))

	return(table_counts)
}

read_table_annotations <- function(filename){
	table_annotations <- read.csv(filename, sep = "\t")

	row.names(table_annotations) <- table_annotations$locusTag
	
	return(table_annotations)
}

read_table_design <- function(filename) {
	table_design <- read.csv(filename, sep = "\t")
	row.names(table_design) <- table_design$sampleId
	
	return(table_design)
}
get_sample_labels <- function(labels){
	#split_labels <- strsplit(labels, "_")
	split_labels <-c()
	
	for (label in labels){
		split_label <- strsplit(label, "_")
		value = split_label <- split_label[[1]][[1]]
		
		split_labels <- c(split_labels, value)
	}
	split_labels <- unique(split_labels)
	return(split_labels)
}
save_deseq_results <- function(results, filename) {
	write.table(as.data.frame(results), filename, sep = "\t", row.names = TRUE)
}

# Includes writing the pairwise as .tsv files, and passing along a master list
run_deseq_pairwise <- function(deseq_object, output_folder, labels) {
  strains <- unique(labels)
  perm <- length(strains)*(length(strains)-1)
  print(strains)
  frames <- vector("list", perm)
  lefts <- vector("list", perm)
  rights <- vector("list", perm)
  count <- 1
	for (left in strains) {
		for (right in strains){
			if (left!=right){
				deseq_results_strain <- results(deseq_object, contrast=c("strain",left, right))
				combo <- paste(left, right, sep="_")
				
				filename = paste(output_folder,"/", left,".",right, ".differentialexpression.tsv", sep = "")
				# print(filename)
				save_deseq_results(deseq_results_strain, filename)
        frames[[count]] <- deseq_results_strain
        lefts[[count]] <- left
        rights[[count]] <- right
        count <- count+1
        #Truthfully probably only needed a list with named elements... look into this to make life easier?
			}
		}
	}
  allByAll <- cbind(frames, lefts, rights)
  colnames(allByAll)<- c("des_result", "strainCompared", "comparedTo")
  return(allByAll)
}
#Creates a deseq_object with the reference as each unique strain
run_deseq_byStrain <- function(deseq_matrix){
  strains <- unique(deseq_matrix$strain)
  all_dobjects <- vector("list", length(strains))
  count <- 1
  for (comparison in strains){
    deseq_matrix$strain <- relevel(deseq_matrix$strain, ref=comparison)
    all_dobjects[count] <- DESeq(deseq_matrix)
    count <- count+1
  }
  names(all_dobjects) <- paste("Ref", strains, sep="")
  return(all_dobjects)
}
run_deseq <- function(table_counts, table_design, table_annotations, output_folder){
	print("Generating the Deseq matrix object...")
	deseq_matrix <- DESeqDataSetFromMatrix(
		countData = table_counts,
		colData = table_design,
		design = ~ strain # Formula summarizing the experimental design
	)
	
	# Make sure the WT is considered the reference sample.
	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "WT")
	print("Generating the deseq object...")
	deseq_object <- DESeq(deseq_matrix)

	print("Running pairwise workflow...")
	deseq_results <- run_deseq_pairwise(deseq_object, output_folder, table_design$strain)
	print("Generating results...")

	# filename_deseq_results <- paste(output_folder, "unknownwithWTascomparison.tsv", sep  = "/")
	# print("Saving results...")
	# save_deseq_results(deseq_results, filename_deseq_results) # deseq_results is now very large
	# Write the deseq results to a table.
	
	
	#Now generatre multiple deseq objects, each with a different reference
	all_deseq_objects <- run_deseq_byStrain(deseq_matrix)
	
	newlist <- list(all_deseq_objects, deseq_results)
	return(newlist)
}



generate_figure_clusters <- function(deseq_object, filename) {
	deseq_object <- estimateSizeFactors (deseq_object)
	deseq_rlog <- rlog ( deseq_object , blind = TRUE )
	rlog.norm.counts <- assay (deseq_rlog)
	distance.m_rlog <- as.dist(1 - cor( rlog.norm.counts, method = "pearson")) #This isn't saving or generating an object
}

generate_figure_pca <- function(deseq_object) {
	deseq_object <- estimateSizeFactors (deseq_object)
	deseq_rlog <- rlog ( deseq_object , blind = TRUE )
	P <- plotPCA(deseq_rlog, intgroup = "strain")
	P + theme_bw() + ggtitle ( " Rlog transformed counts " )
}

# generate_report <- function(deseq_object, table_counts, table_design){
# 	library("ReportingTools")
# 	
# 	desReport <- HTMLReport(
# 		shortName = 'RNAseq_analysis_with_DESeq',
# 		title = 'RNA-seq analysis of differential expression using DESeq, control is WT',
# 		reportDirectory = "./reports"
# 	)
# 	publish(
# 		deseq_object,
# 		desReport,
# 		countTable=table_counts, 
# 		pvalueCutoff=0.05,
# 		factor = table_design$strain,
# 		#conditions = table_annotations$strain,
# 		#conditions=conditions,
# 		#.modifyDF=makeDESeqDF,
# 		expName="deseq",
# 		reportDir="./reports"
# 	)
# 	finish(desReport)
# }

run_pca_interactive <- function(deseq_object, table_design){
	library(pcaExplorer)
	pcaExplorer(dds = deseq_object, coldata = table_design)
}

main <- function(){
	project_folder <- here()
	data_folder <- paste(project_folder, "data", sep = "/")
	results_folder <- paste(project_folder, "results", sep = "/")
	
	
	filename_counts <- paste(data_folder, "abundance.all.matrix.tsv", sep = "/")
	filename_design <- paste(data_folder, "abundance.all.design.tsv", sep = "/")
	filename_annotations <- paste(data_folder, "annotations.tsv", sep = "/")
	
	
		# Specify where the deseq tables should be saved
	folder_deseq <- paste(results_folder, "DESeq2_results_Expanded", sep = "/")
	
	folder_figures <- paste(project_folder, "figures", sep = "/")
	filename_clusters <- paste(folder_figures, "clusters.histogram.png", sep = "/")
	
	dir.create(file.path(folder_deseq), showWarnings = FALSE)
	
	print("Reading in the source tables...")
	# Matrix of read counts
	print("Reading in the count matrix...")
	table_counts <- read_table_counts(filename_counts)
	print("Reading the design table...")
	table_design <- read_table_design(filename_design)
	print("Reading the annotations table...")
	table_annotations <- read_table_annotations(filename_annotations)
	
	print("Running deseq...")
	deseq_tables <- run_deseq(table_counts, table_design, table_annotations, folder_deseq)
	deseq_object <- deseq_tables[[1]]
	deseq_results <- deseq_tables[[2]]
	print("Generating the cluster image...")
	# generate_figure_clusters(deseq_object, filename_clusters) #I don't think these are saved, an will no longer work now that deseq_object is multiples
	
  #Will need to handle LFCshrink and pca in the main analysis part.
	
	newlist <- list(table_counts, table_design, table_annotations, deseq_object, deseq_results)
	return(newlist)
}
							
results <- main()
table_annotations <- results[[3]]
deseq_object <- results[[4]]
table_counts <- results[[1]]
table_design <- results[[2]]
deseq_results <- results[[5]]
newNames <- paste(deseq_results[ , "strainCompared"], deseq_results[ , "comparedTo"], sep="_to_")
dimnames(deseq_results)[[1]] <- newNames
