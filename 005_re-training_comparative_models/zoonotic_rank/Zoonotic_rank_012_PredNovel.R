## parameter

options(error = recover)
options(expressions=50000)

n_models <- 100

format <- "genbank"

output_dir <- commandArgs(trailingOnly=TRUE)[1]
training_data <- commandArgs(trailingOnly=TRUE)[2]

seq_f <- commandArgs(trailingOnly=TRUE)[3]
metadata_f <- commandArgs(trailingOnly=TRUE)[4]

output_f <- "test"

## train modelのあるフォルダ
WORK_DIR <- paste(output_dir, "zoonotic_rank", sep="/")

##
# install.packages("ggplot2", dependencies = TRUE)

suppressPackageStartupMessages({
  library(rprojroot)
  library(seqinr)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(readr)
	library(ggplot2)
  library(xgboost)
  library(caret)
	library(ModelMetrics)
  library(apcluster)
  library(ape)
	library(EnvStats)
})

##
RUN_ID <- paste0("RunID")  	# RunID of trained models to use (only genome feature and mimicry-based runs implemented)
USE_DERIVED_FEATURES <- TRUE  					# Whether derived (mimicry) features should be calculated
MODEL_TYPE <- "xgbTree"       					# Type of model used in this run

##
ROOT_DIR <- "zoonotic_rank"
genome_feature_script <- file.path(ROOT_DIR, "Utils", "GenomeFeatures.py")

source(file.path(ROOT_DIR, "Utils", "cds_parser.R"))
source(file.path(ROOT_DIR, "Utils", "xgboost_utils.R"))
source(file.path(ROOT_DIR, "Utils", "plot_utils.R"))
source(file.path(ROOT_DIR, "Utils", "prediction_utils.R"))
source(file.path(ROOT_DIR, 'Utils', 'calibration_utils.R'))
source(file.path(ROOT_DIR, 'Utils', 'rundata_utils.R'))

if (USE_DERIVED_FEATURES)
	source(file.path(ROOT_DIR, 'Utils', 'derived_genome_feature_utils.R'))


# Trained models:
trainedModels <- paste0(RUN_ID, "_ModelFits.rds") %>%
  file.path(WORK_DIR, "RunData", RUN_ID, .) %>%
  readRDS()

message("\n\nReading models")

## 
if (n_models > length(trainedModels))
  stop(sprintf("n_models is too high: only %d models have been trained in this model set", 
               length(trainedModels)))
 
 ##
 # Predictions from trained models (used to assess accuracy and for calibration)
all_model_preds <- paste0(RUN_ID, "_Predictions.csv") %>% 
	file.path(WORK_DIR, "RunData", RUN_ID, .) %>%
	read.csv()

test_preds <- all_model_preds %>% 
	filter(Dataset == 'test')

calibration_preds <- all_model_preds %>% 
	filter(Dataset == 'calibration')
	
##
# Training data and features
trainingData <- paste0(training_data) %>% 
file.path(WORK_DIR, .) %>%
  read.csv() %>%
  ungroup()

Calculated_dir <- paste0("CalculatedData")

trainingFeatures <- file.path(WORK_DIR, Calculated_dir, "GenomicFeatures-Virus.rds") %>%
	readRDS() %>%
  rename_at(vars(-UniversalName, -Strain), ~ paste0("VirusDirect_", .)) %>%
  filter(paste(UniversalName, Strain) %in% paste(trainingData$UniversalName, trainingData$Strain))

if (USE_DERIVED_FEATURES) {
	# Data needed to calculate derived features:
	humanFeatures <- readRDS(file.path(WORK_DIR, Calculated_dir, 'GenomicFeatures-HumanCombined.rds'))
	
	# Data needed to identify training viruses by feature value:
	virusDistFeatures <- readRDS(file.path(WORK_DIR, Calculated_dir, 'GenomicFeatures-Distances.rds'))
	trainingFeatures <- trainingFeatures %>% 
		left_join(virusDistFeatures, by = c('UniversalName', 'Strain'))
}

print (file.path(WORK_DIR, Calculated_dir, 'GenomicFeatures-HumanCombined.rds'))
print (colnames(humanFeatures))

##
# Metadata for novel viruses:
message("\n\nReading metadata")
print (metadata_f)

expectedCols <- c("Name", "SequenceID", "Start", "Stop")
colTypes <- list(Name = "c", SequenceID = "c", Start = "d", Stop = "d")

if (format == "genbank") {
  expectedCols <- expectedCols[1:2]
  colTypes <- colTypes[1:2]
}

colTypes <- do.call(cols, colTypes)

# Ignore actual header so column names are predictable (but assume there is one):
metaData <- read_csv(metadata_f, col_types = colTypes) # col_names = expectedCols,, skip = 1)

exclude <- NULL

if (!is.null(exclude)) {
	excluded_spp <- str_split(exclude, pattern = ",")[[1]]
	
	if (!all(excluded_spp %in% trainingData$LatestSppName)) {
		not_found <- excluded_spp[!excluded_spp %in% trainingData$LatestSppName]
		
		stop("Some species marked for exclusion were either not in the training data or their names are invalid:\n\t", 
						paste(unique(not_found), collapse = "\n\t"))
	}
}

##
if (format == "fasta") {
    message("\nParsing coding sequences from fasta file\n")

  sequences <- read.fasta(seq_f)

  if (!all(metaData$SequenceID %in% names(sequences)))
    stop(paste("Not all sequence IDs were found in the supplied fasta file. Ensure the second",
               "column of the supplied metadata contains valid IDs (see PredictNovel.R --help)"))

  sequences <- sequences[metaData$SequenceID]
  stopifnot(length(sequences) == nrow(metaData))

  coding <- mapply(extract_cds,
                   sequence = sequences,
                   start_coordinate = metaData$Start,
                   stop_coordinate = metaData$Stop,
  								 MoreArgs = list(allow_complementary = TRUE))
    
  # Output: Each CDS simply gets the same name as the parent sequence (but a different SeqID)
  new_ids <- paste(names(coding), seq_len(length(coding)), sep = "_")
  new_metadata <- metaData %>%
    mutate(SequenceID = new_ids)
  
  names(coding) <- new_ids
  
  temp_coding_file <- tempfile()

  write.fasta(coding, names = names(coding), file.out = temp_coding_file)
}

##
# Always creating a tempfile here, since GenomeFeatures.py expects specific column names
temp_metadata_file <- file.path(WORK_DIR, "_tmp1.csv")

metaData <- metaData %>% 
  distinct(Name, SequenceID)

write.csv(metaData, temp_metadata_file)

print (temp_metadata_file)

##
message("\nCalculating genome features")

temp_out_coding <- file.path(WORK_DIR, "_unknown_coding_temp")
temp_out_entire <- file.path(WORK_DIR, "_unknown_entire_temp")

if (format == "fasta") {
  temp_extra_metadata_file <- tempfile()
  write.csv(new_metadata, temp_extra_metadata_file)
  
  cmd_string_coding <- paste(genome_feature_script, temp_extra_metadata_file, temp_coding_file, 
  													 "--out", temp_out_coding)
  
  cmd_string_entire <- paste(genome_feature_script, temp_metadata_file, seq_f,
  													 "--noncoding", "--out", temp_out_entire)
} else {
  # Rely on GenomeFeatures.py to extract coding sequences from the genbank file
	cmd_string_coding <- paste(genome_feature_script, metadata_f, seq_f,
														 "--extract", "--out", temp_out_coding)
	
	cmd_string_entire <- paste(genome_feature_script, metadata_f, seq_f,
														 "--extract", "--noncoding", "--out", temp_out_entire)
}

# Run these calls:
system2("python3", cmd_string_coding)
system2("python3", cmd_string_entire)

print (temp_out_coding)
print (temp_out_entire)

# Read output
# Currently CPB_machine"s output has some issues:
# (all rows except the header end in a tab, causing data to be shifted relative to column names)
message("\nParsing genome features\n")

features_coding <- read.delim(temp_out_coding, na.strings = c("?", "NA"), row.names = NULL,
															stringsAsFactors = F)
															
colnames(features_coding) <- c(colnames(features_coding)[-1], "XX_REMOVE")
features_coding <- features_coding[, colnames(features_coding) != "XX_REMOVE"]


features_entire <- read.csv(temp_out_entire, stringsAsFactors = F)

# Clean up tempfiles:
# - Capturing return code to preven errors when a file is missing
if (format == "fasta") {
  unlink(temp_coding_file)
  unlink(temp_extra_metadata_file)
}

# Merge feature sets
features_coding <- features_coding %>% 
	select(-TaxID, -Species) %>% 
	rename_at(vars(-SeqName), ~ paste(., 'Coding', sep = '_'))

features_entire <- features_entire %>% 
	rename_at(vars(-SeqName), ~ paste(., 'EntireSeq', sep = '_'))

all_direct_features <- features_coding %>% 
	full_join(features_entire, by = 'SeqName')
	
# Rename direct features to their final names: 
genomeFeatures <- all_direct_features %>% 
	rename_at(vars(-SeqName), ~ paste("VirusDirect", ., sep = '_')) %>% 
	rename(Name = SeqName)
	
# Calculate derived genome features ('distance' from human genes) if needed
if (USE_DERIVED_FEATURES) {
	feature_names <- colnames(humanFeatures)
	feature_names <- feature_names[! feature_names %in% c('GeneID', 'meanCPM', 'TranscriptID', 'GeneSet')]
	
	all_direct_features <- all_direct_features %>% 
		mutate(UniversalName = SeqName,
					 Strain = NA_character_)  # These columns expected by get_feature_dists()
	
	isg_features <- humanFeatures %>% 
		filter(GeneSet == 'ISG')
	
	isg_features <- get_feature_dists(virusFeatures = all_direct_features, 
																		geneFeatures = isg_features,
																		setprefix = 'ISG',
																		featureColNames = feature_names)
	
	housekeeping_features <- humanFeatures %>% 
		filter(GeneSet == 'Housekeeping')
	
	housekeeping_features <- get_feature_dists(virusFeatures = all_direct_features, 
																						 geneFeatures = housekeeping_features,
																						 setprefix = 'Housekeeping',
																						 featureColNames = feature_names)
	
	remaining_features <- humanFeatures %>% 
		filter(GeneSet == 'Remaining')
	
	remaining_features <- get_feature_dists(virusFeatures = all_direct_features, 
																					geneFeatures = remaining_features,
																					setprefix = 'Remaining',
																					featureColNames = feature_names)
	
	# Join:
	all_indirect_features <- isg_features %>% 
		full_join(housekeeping_features, by = c('UniversalName', 'Strain')) %>% 
		full_join(remaining_features, by = c('UniversalName', 'Strain')) %>% 
		select(-Strain)
	
	genomeFeatures <- genomeFeatures %>% 
		left_join(all_indirect_features, by = c('Name' = 'UniversalName'))
}


##
if (!is.null(exclude)) {
	virus_names <- trainingData %>% 
		distinct(UniversalName, LatestSppName)
	
	test_preds <- test_preds %>% 
		left_join(virus_names, by = "UniversalName")
	
	valid_models <- test_preds %>% 
		group_by(Iteration) %>% 
		summarise(keep = all(excluded_spp %in% LatestSppName)) %>% 
		filter(keep) %>% 
		pull(Iteration)
	
	if (length(valid_models) < n_models)
		stop("Not enough models remain after excluding the requested viruses. ", 
				 "Adjust either 'n_models' or the exclusion list.")
	
	cat(length(valid_models), "models remain after excluding requested viruses\n")
	
	test_preds <- test_preds %>% 
		filter(Iteration %in% valid_models) %>% 
		filter(!LatestSppName %in% excluded_spp) %>% 
		select(-LatestSppName)
}

##
message("\nCollecting top models\n")

trainedAccuracy <- test_preds %>% 
	mutate(InfectsHumans = InfectsHumans == 'True') %>% 
	group_by(Iteration) %>% 
	summarise(AUC = ModelMetrics::auc(actual = InfectsHumans, predicted = RawScore)) %>% 
	ungroup() %>% 
  mutate(Rank = rank(-AUC, ties.method = "random")) 

topModels <- trainedAccuracy %>% 
  filter(Rank <= n_models) %>% 
  pull(Iteration)

print (topModels)

bagModels <- trainedModels[topModels]

# Predict:
message("\nPrediction\n")

predictionData <- lapply(bagModels, prepare_prediction_data, newdata = genomeFeatures)

allPredictions <- mapply(get_prediction, fit = bagModels, newdata = predictionData,
                         MoreArgs = list(virusnames = genomeFeatures$Name,
                         								 model_type = MODEL_TYPE),
                         SIMPLIFY = FALSE, USE.NAMES = TRUE) %>%
  bind_rows()
  
 ##
 allPredictions <- allPredictions %>% 
	select(.data$Iteration, .data$Name, RawScore = .data$True) %>% 
	group_by(.data$Iteration) %>% 
	group_modify(~ calibrate_preds(.x, calibration_preds = calibration_preds))

##
outDir <- dirname(output_f)

if (!dir.exists(outDir))
  dir.create(outDir, recursive = TRUE)

## Predictions ------------------------------------------------------------------------------------
write_excel_csv(allPredictions, path = sprintf("%s.predictions.csv", output_f))
