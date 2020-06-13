#Requires DESeq2, here, and dplyr libraries loaded.

read_table_counts <- function(filename){
  table_counts <- read.csv(filename, sep = "\t", row.names = 'locusTag')
  # Need to convert the table counts to integers because they may be float with zeros after the decimal point.	
  #table_counts[] <- lapply(table_counts, function(x) as.numeric(x))
  
  return(table_counts)
}
filtering <- function(df, exclusions){
  # For filtering TilS by chromosome/plasmid
  parts <- list(  chr3=paste0("Bcen2424_", seq(5867,6788)),
                  plmd=paste0("Bcen2424_", seq(6789,6947)),
                  both=paste0("Bcen2424_", seq(5867,6947))
  )
  
  table_counts <- df[ !(row.names(df) %in% parts[[exclusions]]), ]
  
  return(table_counts)
}


# read_table_annotations <- function(filename){
#   table_annotations <- read.csv(filename, sep = "\t")
#   
#   row.names(table_annotations) <- table_annotations$locusTag
#   
#   return(table_annotations)
# }

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
        
        filename = paste(output_folder,"/", left,".to.",right, ".tsv", sep = "")
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
  newNames <- paste(lefts, rights, sep="_to_")
  row.names(allByAll) <- newNames
  return(allByAll)
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
  
  deseq_results <- run_deseq_pairwise(deseq_object, here(output_folder), table_design$strain)
  # filename_deseq_results <- paste(output_folder, "unknownwithWTascomparison.tsv", sep  = "/")
  # print("Saving results...")
  # save_deseq_results(deseq_results, filename_deseq_results) # deseq_results is now very large
  
  newlist <- list(deseq_object, deseq_results)
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


des_workflow <- function(resfolder, count_matrix, exclusions, table_design, table_annotations){
  # Expects datafolder and resfolder pathed from the project folder, and matrix and design tables as data frames,
  # annotationsdf in R data.frame/matrix format,
  # and the exclusions as a choice from chr3, plmd, or both.
  
  project_folder <- here()
  
  # filename_counts <- here(datafolder,matrixtsv)
  # filename_design <- here(datafolder,designtsv)
  # filename_annotations <- here(datafolder, annotationstsv)
  
  # folder_figures <- paste(project_folder, "figures", sep = "/")
  # filename_clusters <- paste(folder_figures, "clusters.histogram.png", sep = "/")
  
  dir.create(here(resfolder), showWarnings = FALSE) 
  
  print("Reading in the source tables...")
  # Matrix of read counts
  # print("Reading in the count matrix...")
  # table_counts <- read_table_counts(filename_counts)
  print("Filtering count matrix...")
  table_counts <- filtering(count_matrix, exclusions) #Toggle this to filter out chr3
  
  # print("Reading the design table...")
  # table_design <- read_table_design(filename_design)
  # print("Reading the annotations table...")
  # table_annotations <- read_table_annotations(filename_annotations)
  
  print("Running deseq...")
  deseq_tables <- run_deseq(table_counts, table_design, table_annotations, resfolder)
  deseq_object <- deseq_tables[[1]]
  deseq_results <- deseq_tables[[2]]
  # print("Generating the cluster image...")
  # generate_figure_clusters(deseq_object, filename_clusters) #I don't think these are saved, an will no longer work now that deseq_object is multiples
  
  #Will need to handle LFCshrink and pca in the main analysis part.
  
  newlist <- list(deseq_object, deseq_results)
  return(newlist)  # Results are a list of the deseq object, results list
}

parse_GB_annotation <- function(filename){
  # Expects the file location of the GenBank _feature_table.txt as a string.
  # Returns the cleaned feature table, original feature table, and indeces of non-unique entries  
  
  # Read in the data frame and reformat
  GB_frame <- read.csv(filename, sep = "\t", stringsAsFactors = FALSE, row.names = NULL)
  names(GB_frame)[names(GB_frame) == "X..feature"] <- "feature"
  names(GB_frame)[names(GB_frame) == "locus_tag"] <- "locusTag"
  GB_frame <- GB_frame %>% dplyr::select("locusTag", everything())
  
  # Initialize lsits and tables
  u_features <- dplyr::slice(GB_frame, 0)
  u_list <- vector("character") # List for storing unique locus_tag
  not_unique <- vector("integer") # List of indeces of main table of repeated locus_tags
  count <- 1 #Assumes at least 1 unique feature
  # not_count <- 1 #Assumes at least 1 non-unique feature
  
  # For each feature listed, compares locus_tag with growing list of already seen locus tags,
  # ensures has a name, and copies unique features, and notes where duplicates exist
  for (i in 1:nrow(GB_frame)){
    if ( !(GB_frame[i, "locusTag"] %in% u_list) && (
      (!(GB_frame[i, "name"] == "") | (GB_frame[i, "attributes"] =="pseudo")) )
    ){
      u_features[count, ] <- GB_frame[i, ]
      u_list[count] <- GB_frame[i, "locusTag"]
      count <- count+1
      # } else {
      #   not_unique[not_count] <- i
      #   not_count <- not_count+1
    }
  }
  
  # Format the data table
  rownames(u_features) <- u_features[ ,1]
  for (i in rownames(u_features)) {
    if (u_features[i, "class"] == "") {
      u_features[i, "class"] <- u_features[i, "feature"]
    }
  }
  
  # Save tab-delineated in same folder
  newname <- regmatches(filename, regexpr(".*/", filename))
  newfile <- paste0(newname, "GB.annotations.tsv")
  write.table(u_features, file = newfile, row.names = FALSE, sep = "\t", na = "")
  
  return (u_features)
}

parse_RS_annotation <- function(filename){
  # Expects the file location of the RefSeq _feature_table.txt as a string.
  # Returns the cleaned feature table and original feature table.
  
  # Read in the data frame and reformat
  RS_frame <- read.csv(filename, sep = "\t", stringsAsFactors = FALSE, row.names = NULL)
  names(RS_frame)[names(RS_frame) == "X..feature"] <- "feature"
  names(RS_frame)[names(RS_frame) == "locus_tag"] <- "locusTag"
  RS_frame <- RS_frame %>% dplyr::select("locusTag", everything())
  
  # The RefSeq table had evenly doubled entries. Checking to see which has more info.
  splits <- split(RS_frame, 1:2)
  
  if (any(splits[[1]]$name !="")){
    ind <- 1
    otherind <-2
  } else {
    ind <- 2
    otherind <-1
  }
  
  # Format the data table
  rownames(splits[[ind]]) <- splits[[ind]][ ,1]
  rownames(splits[[otherind]]) <- splits[[ind]][ ,1]
  for (i in rownames(splits[[ind]])) {
    if(splits[[ind]][i, "feature"] == "tRNA"){
      splits[[ind]][i, "attributes"] <- paste0(splits[[ind]][i, "attributes"], "; ", splits[[otherind]][i,"attributes"])
    } else {  
      splits[[ind]][i, "attributes"] <- splits[[otherind]][i, "attributes"]
    }
  }
  anno_table <- splits[[ind]]
  
  gbLocusTag <- vector("character", nrow(anno_table))
  gbLocusTag <- lapply( row.names(anno_table), function(listitem){
    query <- anno_table[listitem, "attributes"]
    if ( query == "") ntag <- anno_table[listitem, "locusTag"]
    else {
      temp <- regmatches(query, regexpr("Bcen2424_\\d{4}", query))
      if ( identical(temp, character(0)) ) ntag <- anno_table[listitem, "locusTag"]
      else ntag <- temp
    }
    return(ntag)
  })
  
  gbLocusTag <- unlist(gbLocusTag)
  names(gbLocusTag) <- row.names(anno_table)
  anno_table <- cbind(anno_table, gbLocusTag)
  names(anno_table)[names(anno_table) == "locusTag"] <- "RSlocusTag"
  names(anno_table)[names(anno_table) == "gbLocusTag"] <- "locusTag"
  anno_table <- anno_table %>% dplyr::select("locusTag", everything())
  row.names(anno_table) <- anno_table[ ,1]
  
  # Save tab-delineated in same folder
  newname <- regmatches(filename, regexpr(".*/", filename))
  newfile <- paste0(newname, "RS.annotations.tsv")
  write.table(anno_table, file = newfile, row.names = FALSE, sep = "\t", na = "")
  return (anno_table)
}

clustering_plot <- function(sampledists, rld, colors, title, filepathname){
# Requires the sample distances, transformed data, color scheme, title of plot, and path/name of file to save including .pdf
  sdMatrix <- as.matrix(sampledists)
  rownames(sdMatrix) <- paste(rld$strain, rld$replicate, sep="-")
  colnames(sdMatrix) <- NULL
  pdf(file = filepathname) # Uncomment this and the dev.off to save individual files
  pheat <- pheatmap(sdMatrix,
                         clustering_distance_rows=sampledists,
                         clustering_distance_cols=sampledists,
                         col=colors, main = title)
  dev.off()
  return(pheat)
}
