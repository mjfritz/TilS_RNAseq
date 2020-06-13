library("DESeq2")

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

# run_deseq_pairwise <- function(deseq_object, output_folder, labels) {
# 	for (left in labels) {
# 		for (right in labels){
# 			if (left!=right){
# 				deseq_results_strain <- results(deseq_object, contrast=c("strain",left, right))
# 				
# 				filename = paste(output_folder,"/", left,".",right, ".differentialexpression.tsv", sep = "")
# 				print(filename)
# 				save_deseq_results(deseq_results_strain, filename)
# 			}
# 		}
# 	}
# }

run_deseq <- function(table_counts, table_design, table_annotations, output_folder){
	print("\tGenerating the Deseq matrix object...")
	deseq_matrix <- DESeqDataSetFromMatrix(
		countData = table_counts,
		colData = table_design,
		design = ~ strain # Formula summarizing the experimental design
	)
	
	# for(str in deseq_matrix$strain){
	#   cat("Generating the ", str, " deseq object and results.")
	#   strObject <- paste("deseq_ojbect",str,sep="_")
	#   deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = str)
	#   assign(strObject, DESeq(deseq_matrix))
	#   get(strObject)
	#   Paste("deseq_results",str,sep="_") <- results(get(strObject)) # This line is now problematic
	# }


	# Running deseq with each strain as the control, first WT, then beyond.
	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "WT")
	#deseq_matrix$ppc <- relevel(deseq_matrix$group, ref = 'group')
	print("Generating the WT deseq object...")
	deseq_object_WT <- DESeq(deseq_matrix)
	deseq_results_WT <- results(deseq_object_WT) #These are only generating tRNA against the strain

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "A244T")
	deseq_object_A244T <- DESeq(deseq_matrix)
	deseq_results_A244T <- results(deseq_object_A244T)

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "N274Y")
	deseq_object_N274Y <- DESeq(deseq_matrix)
	deseq_results_N274Y <- results(deseq_object_N274Y)

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "N445K")
	deseq_object_N445K <- DESeq(deseq_matrix)
	deseq_results_N445K <- results(deseq_object_N445K)

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "P421L")
	deseq_object_P421L <- DESeq(deseq_matrix)
	deseq_results_P421L <- results(deseq_object_P421L)

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "S3C")
	deseq_object_S3C <- DESeq(deseq_matrix)
	deseq_results_S3C <- results(deseq_object_S3C)

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "T41T")
	deseq_object_T41T <- DESeq(deseq_matrix)
	deseq_results_T41T <- results(deseq_object_T41T)

	deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "tRNA")
	deseq_object_tRNA <- DESeq(deseq_matrix)
	deseq_results_tRNA <- results(deseq_object_tRNA)

	#strain_order <- c("WT","A244T","N274Y","N445K", "P421L", "S3C", "T41T", "tRNA")
	#strain_order <- get_sample_labels(table_design$strain)
	#print("Running pairwise workflow...")
	#deseq_results <- run_deseq_pairwise(deseq_object, output_folder, table_design$strain)
	print("Generating results...")

	filename_deseq_results <- paste(output_folder, "WT.differentialexpression.tsv", sep  = "/") 
	save_deseq_results(deseq_results_WT, filename_deseq_results) # Use of naming variable can be ... expanded?
	filename_deseq_results <- paste(output_folder, "A244T.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_A244T, filename_deseq_results)
	filename_deseq_results <- paste(output_folder, "N274Y.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_N274Y, filename_deseq_results)
	filename_deseq_results <- paste(output_folder, "N445K.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_N445K, filename_deseq_results)
	filename_deseq_results <- paste(output_folder, "P421L.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_P421L, filename_deseq_results)
	filename_deseq_results <- paste(output_folder, "S3C.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_S3C, filename_deseq_results)
	filename_deseq_results <- paste(output_folder, "T41T.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_T41T, filename_deseq_results)
	filename_deseq_results <- paste(output_folder, "tRNA.differentialexpression.tsv", sep  = "/")
	save_deseq_results(deseq_results_tRNA, filename_deseq_results)
	
	
	# Write the deseq results to a table.
	newlist <- list(deseq_object_WT, deseq_results_WT,
	                deseq_object_A244T, deseq_results_A244T,
	                deseq_object_N274Y, deseq_results_N274Y,
	                deseq_object_N445K, deseq_results_N445K,
	                deseq_object_P421L, deseq_results_P421L,
	                deseq_object_S3C, deseq_results_S3C,
	                deseq_object_T41T, deseq_results_T41T,
	                deseq_object_tRNA, deseq_results_tRNA) #I think this passes everything on as a single object
}



generate_figure_clusters <- function(deseq_object, filename) {
	deseq_object <- estimateSizeFactors (deseq_object)
	deseq_rlog <- rlog ( deseq_object , blind = TRUE )
	rlog.norm.counts <- assay (deseq_rlog)
	distance.m_rlog <- as.dist(1 - cor( rlog.norm.counts, method = "pearson"))
	#do these actually get saved with the filename variable?
}

generate_figure_pca <- function(deseq_object) {
	deseq_object <- estimateSizeFactors (deseq_object)
	deseq_rlog <- rlog ( deseq_object , blind = TRUE )
	P <- plotPCA(deseq_rlog, intgroup = "strain")
	P + theme_bw() + ggtitle ( " Rlog transformed counts " )
} #Not sure if this gets used?

generate_report <- function(deseq_object, table_counts, table_design){
	library("ReportingTools")
	
	desReport <- HTMLReport(
		shortName = 'RNAseq_analysis_with_DESeq',
		title = 'RNA-seq analysis of differential expression using DESeq',
		reportDirectory = "./reports"
	)
	publish(
		deseq_object,
		desReport,
		countTable=table_counts, 
		pvalueCutoff=0.05,
		factor = table_design$strain,
		#conditions = table_annotations$strain,
		#conditions=conditions,
		#.modifyDF=makeDESeqDF,
		expName="deseq",
		reportDir="./reports"
	)
	finish(desReport)
}

run_pca_interactive <- function(deseq_object, table_design){
	library(pcaExplorer)
	pcaExplorer(dds = deseq_object, coldata = table_design)
}

main <- function(){
	project_folder <- here()
	#project_folder <- "/home/cld100/storage/projects/eisha_rna/"
	data_folder <- paste(project_folder, "data", sep = "/")
	results_folder <- paste(project_folder, "results", sep = "/")
	
	
	filename_counts <- paste(data_folder, "abundance.all.matrix.tsv", sep = "/")
	filename_design <- paste(data_folder, "abundance.all.design.tsv", sep = "/")
	filename_annotations <- paste(data_folder, "annotations.tsv", sep = "/")
	
	
	# Specify where the deseq tables should be saved
	folder_deseq <- paste(results_folder, "DESeq2_results_pairwise", sep = "/")
	
	folder_figures <- paste(project_folder, "figures", sep = "/")
	filename_clusters_WT <- paste(folder_figures, "WT.clusters.histogram.png", sep = "/")
	filename_clusters_A244T <- paste(folder_figures, "A244T.clusters.histogram.png", sep = "/")
	filename_clusters_N274Y <- paste(folder_figures, "N274Y.clusters.histogram.png", sep = "/")
	filename_clusters_N445K <- paste(folder_figures, "N445K.clusters.histogram.png", sep = "/")
	filename_clusters_P421L <- paste(folder_figures, "P421L.clusters.histogram.png", sep = "/")
	filename_clusters_S3C <- paste(folder_figures,"S3C.clustes.histogram.png", sep ="/")
	filename_clusters_T41T <- paste(folder_figures, "T41T.clusters.histogram.png", sep = "/")
	filename_clusters_tRNA <- paste(folder_figures, "tRNA.clusters.histogram.png", sep = "/")

	
	
	dir.create(file.path(folder_deseq), showWarnings = FALSE)

	
	print("Reading in the source tables...")
	# Matrix of read counts
	print("Reading in the count matrix...")
	table_counts <- read_table_counts(filename_counts)
	print("Reading the design table...")
	table_design <- read_table_design(filename_design)
	print("Reading the annotations table...")
	table_annotations <- read_table_annotations(filename_annotations)
	
	#Should have... created matrices in main, and called run_deseq on each to only generate two tables at a time
	print("Running deseq...")
	deseq_tables <- run_deseq(table_counts, table_design, table_annotations, folder_deseq) # Returning list of al of the tables
	deseq_object_WT <- deseq_tables[[1]]
	deseq_results_WT <- deseq_tables[[2]]
	deseq_object_A244T <- deseq_tables[[3]]
	deseq_results_A244T <- deseq_tables[[4]]
	deseq_object_N274Y <- deseq_tables[[5]]
	deseq_results_N274Y <- deseq_tables[[6]]
	deseq_object_N445K <- deseq_tables[[7]]
	deseq_results_N445K <- deseq_tables[[8]]
	deseq_object_P421L <- deseq_tables[[9]]
	deseq_results_P421L <- deseq_tables[[10]]
	deseq_object_S3C <- deseq_tables[[11]]
	deseq_results_S3C <- deseq_tables[[12]]
	deseq_object_T41T <- deseq_tables[[13]]
	deseq_results_T41T <- deseq_tables[[14]]
	deseq_object_tRNA <- deseq_tables[[15]]
	deseq_results_tRNA <- deseq_tables[[16]]
	
	print("Generating the cluster image...")
	generate_figure_clusters(deseq_object_WT, filename_clusters_WT)
	generate_figure_clusters(deseq_object_A244T, filename_clusters_A244T)
	generate_figure_clusters(deseq_object_N274Y, filename_clusters_N274Y)
	generate_figure_clusters(deseq_object_N445K, filename_clusters_N445K)
	generate_figure_clusters(deseq_object_P421L, filename_clusters_P421L)
	generate_figure_clusters(deseq_object_S3C, filename_clusters_S3C)
	generate_figure_clusters(deseq_object_T41T, filename_clusters_T41T)
	generate_figure_clusters(deseq_object_tRNA, filename_clusters_tRNA)
	
	
	#resLFC <- lfcShrink(deseq_object, coef="strain_A244T_vs_WT", type="apeglm") #Need 7 x 7 of these?
	#print("Generating report...")
	#generate_report(deseq_object, table_counts, table_design) #Need 8 of these
	#print("Generating PCA figure...")
	#generate_figure_pca(deseq_object) # Need 8? of these
	
	newlist <- list(table_counts, table_design, table_annotations,
	                deseq_object_WT, deseq_results_WT, #4 and 5
	                deseq_object_A244T, deseq_results_A244T,
	                deseq_object_N274Y, deseq_results_N274Y,
	                deseq_object_N445K, deseq_results_N445K,
	                deseq_object_P421L, deseq_results_P421L,
	                deseq_object_S3C, deseq_results_S3C,#14 and 15
	                deseq_object_T41T, deseq_results_T41T,
	                deseq_object_tRNA, deseq_results_tRNA)
	return(newlist)
}
							
results <- main()
table_annotations <- results[[3]]
table_counts <- results[[1]]
table_design <- results[[2]]
deseq_object_WT <- results[[4]]
deseq_results_WT <- results[[5]]
deseq_object_A244T <- results[[6]]
deseq_results_A244T <- results[[7]]
deseq_object_N274Y <- results[[8]]
deseq_results_N274Y <- results[[9]]
deseq_object_N445K <- results[[10]]
deseq_results_N445K <- results[[11]]
deseq_object_P421L <- results[[12]]
deseq_results_P421L <- results[[13]]
deseq_object_S3C <- results[[14]]
deseq_results_S3C <- results[[15]]
deseq_object_T41T <- results[[16]]
deseq_results_T41T <- results[[17]]
deseq_object_tRNA <- results[[18]]
deseq_results_tRNA <- results[[19]]

# All of these strain specific manipulations could be easier with vectors, I think. But I'm not sure how to do it yet.