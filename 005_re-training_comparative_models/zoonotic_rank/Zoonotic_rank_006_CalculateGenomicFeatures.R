library(tidyr)
library(dplyr)
library(seqinr)
library(EnvStats)
library(matrixStats)

source(file.path('zoonotic_rank', 'Utils', 'derived_genome_feature_utils.R'))

CPM_CUTOFF <- 1 # mean CPM values above this are considered as evidence that the gene is expressed

## input
input_f <- commandArgs(trailingOnly=TRUE)[1]
fasta_f <- commandArgs(trailingOnly=TRUE)[2]
seq_dir <- commandArgs(trailingOnly=TRUE)[3]
output_dir <- commandArgs(trailingOnly=TRUE)[4]

work_dir <- paste(output_dir, "zoonotic_rank", sep="/")

VIRUS_DATA <- paste(input_f, sep="/")
VIRUS_SEQS_GENBANK <- paste(work_dir, "Sequences", sep="/")
VIRUS_SEQS_FASTA <- paste(work_dir, "Sequences", "curated.fasta", sep="/")

calculated_key <- paste0("CalculatedData")
calculated_dir <- paste(work_dir, calculated_key, sep="/")

dir.create(calculated_dir)

## input (fixed)
ISG_IDENTITY <- 'zoonotic_rank/InternalData/Shaw2017_raw/ISG_PublishedData_Web.csv'
EXPRESSION_DATA <- 'zoonotic_rank/InternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv'  # TODO: should be using TPM instead of these CPM values

HOUSEKEEPING_IDENTITY <- 'zoonotic_rank/CalculatedData/HumanGeneSets/HousekeepingGeneIDs.csv'

TRANSCRIPT_DATA <- 'zoonotic_rank/CalculatedData/HumanGeneSets/TranscriptData.csv'
TRANSCRIPT_SEQS <- 'zoonotic_rank/CalculatedData/HumanGeneSets/TranscriptSequences.fasta'

##
DINUCLEOTIDE <- '^[AUGC]p[ATUGC]'
BRIDGE_DINUCLEOTIDE <- '^br[AUGC]p[AUGC]'
NON_BRIDGE_DINUCLEOTIDE <- 'NonBr[AUGC]p[AUGC]'
AA_BIAS <- '^[A-Z].Bias'
CODON_BIAS <- '[ATGC]{3}.Bias'

##
NUCLEOTIDE_CONTENT <- '^[ATGCN]_'
AT_GC_CONTENT <- '^[ATGC]{2}_'
CODON_PAIR_BIAS <- '[ATGC]{3}.[A-Z]..[ATGC]{3}.[A-Z].'
CPS <- 'Cps'  # TODO: Not sure what these two columns are

ALL_FEATURES <- paste(NUCLEOTIDE_CONTENT, AT_GC_CONTENT, DINUCLEOTIDE, BRIDGE_DINUCLEOTIDE,
                      NON_BRIDGE_DINUCLEOTIDE, AA_BIAS, CODON_BIAS, CODON_PAIR_BIAS, CPS,
                      sep = '|')

USED_FEATURES <- paste(DINUCLEOTIDE, BRIDGE_DINUCLEOTIDE, NON_BRIDGE_DINUCLEOTIDE,
                       AA_BIAS, CODON_BIAS, sep = '|')
                       
##
calculate_genomic <- function(metaData, sequenceLocation, extract, noncoding) {
    ## Write data and call python script
    temp_f <- paste(work_dir, "_temp.csv", sep="/")
    write.csv(metaData, temp_f)
    cmdString <- paste('zoonotic_rank/Utils/GenomeFeatures.py', temp_f,
                       sequenceLocation)
    if (extract) cmdString <- paste(cmdString, '--extract')
    if (noncoding) cmdString <- paste(cmdString, '--noncoding')
    
    system2('python3', cmdString) ## ここに問題あり
    
    ## Read in result
    if (!noncoding) {
        # Currently CPB_machine's output has some issues:
        # - all rows except the header end in a tab, causing data to be shifted relative to column names:
        temp_genome_f <- paste(work_dir, "_temp_GenomeFeatures.csv", sep="/")
        result <- read.delim(temp_genome_f, na.strings = c('?', 'NA'),
                             row.names = NULL,
                             stringsAsFactors = F)
        colnames(result) <- c(colnames(result)[-1], 'XX_REMOVE')
        result <- result[, colnames(result) != 'XX_REMOVE']
    } else {
        temp_genome_f <- paste(work_dir, "_temp_GenomeFeatures.csv", sep="/")
        result <- read.csv(temp_genome_f)
    }
    
    ## Clean up temp files and return
#     unlink('result/20230404/Filoviridae/zoonotic_rank/_temp.csv')
#     unlink('result/20230404/Filoviridae/zoonotic_rank/_temp_GenomeFeatures.csv')
    return(result)
}

##
check_missing_seqs <- function(transcriptData, transcriptSeqs) {
    seqData <- names(transcriptSeqs) %>%
    data.frame(SeqName = ., stringsAsFactors = F) %>%
    separate(SeqName, into = c('GeneID', 'TranscriptID'), sep = '_')
    
    if (any(is.na(transcriptData$TranscriptID)))
        warning('Some genes are missing matching transcript data')
    
    if (! all(transcriptData$TranscriptID %in% seqData$TranscriptID) )
        warning('Some transcript sequences are missing')
    
    transcriptData %>%
    filter(! TranscriptID %in% seqData$TranscriptID)
}

##
virusMetaData <- read.csv(VIRUS_DATA)
isgIdentityData <- read.csv(ISG_IDENTITY, stringsAsFactors = FALSE)
ExpressionData <- read.csv(EXPRESSION_DATA, stringsAsFactors=FALSE)
colnames(ExpressionData) <- c('GeneID', colnames(ExpressionData)[-1]) # First column unnamed in file
housekeepingData <- read.csv(HOUSEKEEPING_IDENTITY)
transcriptData <- read.csv(TRANSCRIPT_DATA)
transcriptSequences <- read.fasta(TRANSCRIPT_SEQS)

##
virusMetaData <- virusMetaData %>%
    select(UniversalName, Strain, Accessions) %>%
    mutate(virusID = paste(UniversalName, Strain, sep = '_'))
    
##
featureDat <- virusMetaData %>%
    select(Name = virusID,
           SequenceID = Accessions) %>%
    mutate(SequenceID = strsplit(as.character(SequenceID), split = "; ")) %>%
    unnest(cols = c(SequenceID))

featureDat <- data.frame(lapply(featureDat, as.character), stringsAsFactors = FALSE)

print (head(featureDat))

virusCoding <- calculate_genomic(featureDat, VIRUS_SEQS_GENBANK, extract = TRUE, noncoding=FALSE) %>%
    select(SeqName, matches(USED_FEATURES)) %>%
    select(-X.Bias, -ATG.Bias, -TGG.Bias) %>%   # These have 0 variance (while X bias is not a real feature)
    rename_at(vars(-SeqName), ~ paste(., 'Coding', sep = '_'))

##
temp_prob_f <- paste(work_dir, "_temp_Problems.csv", sep="/")
feature_f <- paste(calculated_dir, "GenomicFeatures_CDSExtractionFailed-Virus.csv", sep="/")
file.rename(temp_prob_f, feature_f)

##
virusEntireSeq <- calculate_genomic(featureDat, VIRUS_SEQS_FASTA, extract=FALSE, noncoding = TRUE) %>%
    select(SeqName, matches(USED_FEATURES)) %>%
    rename_at(vars(-SeqName), ~ paste(., 'EntireSeq', sep = '_'))


virusFeatures <- virusCoding %>%
    full_join(virusEntireSeq, by = 'SeqName') %>%
    left_join(virusMetaData, by = c('SeqName' = 'virusID'))

##
transcriptData <- transcriptData %>%
    mutate(SeqName = paste(GeneID, TranscriptID, sep = '_')) %>%
    filter(SeqName %in% names(transcriptSequences)) %>%
    distinct()

codingFeatures <- transcriptData %>%
    filter(Biotype == 'protein_coding') %>%
    select(SequenceID = SeqName) %>%
    mutate(Name = SequenceID) %>%
    calculate_genomic(TRANSCRIPT_SEQS, extract=FALSE, noncoding = FALSE) %>%
    filter(!is.na(A)) %>%  # See TODO above
    select(SeqName, matches(USED_FEATURES)) %>%
    rename_at(vars(-SeqName), ~ paste(., 'Coding', sep = '_'))

temp_prob_f <- paste(work_dir, "_temp_Problems.csv", sep="/")
feature_f <- paste(calculated_dir, "GenomicFeatures_CDSExtractionFailed-Human.csv", sep="/")
file.rename(temp_prob_f, feature_f)

##
noncodingFeatures <- transcriptData %>%
    select(SequenceID = SeqName) %>%
    mutate(Name = SequenceID) %>%
    calculate_genomic(TRANSCRIPT_SEQS, extract=FALSE, noncoding = TRUE) %>%
    select(SeqName, matches(USED_FEATURES)) %>%
    rename_at(vars(-SeqName), ~ paste(., 'EntireSeq', sep = '_'))

combinedFeatures <- codingFeatures %>%
    full_join(noncodingFeatures, by = 'SeqName') %>%
    separate(SeqName, into = c('GeneID', 'TranscriptID'), sep = '_')
    
##
ExpressionData <- ExpressionData %>%
    gather(-GeneID, key = 'key', value = 'CPM') %>%
    separate(key, into = c('Condition', 'Replicate')) %>%
    group_by(GeneID, Condition) %>%
    summarise(meanCPM = mean(CPM)) %>%
    ungroup() %>%
    filter(meanCPM >= CPM_CUTOFF) %>%
    filter(! startsWith(GeneID, 'CVR'))

MockExpression <- ExpressionData %>%
    filter(Condition == 'C1')

StimulatedExpression <- ExpressionData %>%
    filter(Condition == 'C2')
    
isgIDs <- isgIdentityData %>%
    filter(Species == 'Homo sapiens') %>%
    filter(Expression == 'up_regulated') %>%
    filter(! startsWith(ENSEMBL.ID, 'CVR')) %>%
    .$ENSEMBL.ID

housekeepingIDs <- housekeepingData %>%
    filter(! GeneID %in% isgIDs) %>%   									# Removes 321 housekeeping genes which are interferon-stimulated
    filter(GeneID %in% MockExpression$GeneID) %>%  	# Removes 97 genes not detected in Shaw et al's data (at least not in mock-stimulated data)
    .$GeneID

remainingIDs <- MockExpression %>%
    filter(! GeneID %in% isgIDs) %>%
    filter(! GeneID %in% housekeepingIDs) %>%
    filter(! startsWith(GeneID, 'CVR')) %>%
    .$GeneID


isgFeatures <- StimulatedExpression %>%
    filter(GeneID %in% isgIDs) %>%
    left_join(combinedFeatures, by = 'GeneID') %>%
    ungroup()

housekeepingFeatures <- MockExpression %>%
    filter(GeneID %in% housekeepingIDs) %>%
    left_join(combinedFeatures, by = 'GeneID') %>%
    ungroup()

remainingFeatures <- MockExpression %>%
    filter(GeneID %in% remainingIDs) %>%
    left_join(combinedFeatures, by = 'GeneID') %>%
    ungroup()
    
stopifnot(nrow(isgFeatures) == length(unique(isgFeatures$GeneID)))
stopifnot(nrow(housekeepingFeatures) == length(unique(housekeepingFeatures$GeneID)))
stopifnot(nrow(remainingFeatures) == length(unique(remainingFeatures$GeneID)))

##
featureColNames <- colnames(virusFeatures)[grepl(ALL_FEATURES, colnames(virusFeatures))]

isgDists <- get_feature_dists(virusFeatures,
                              geneFeatures = isgFeatures,
                              setprefix = 'ISG',
                              featureColNames = featureColNames)

housekeepingDists <- get_feature_dists(virusFeatures,
                                       geneFeatures = housekeepingFeatures,
                                       setprefix = 'Housekeeping',
                                       featureColNames = featureColNames)

remainingDists <- get_feature_dists(virusFeatures,
                                    geneFeatures = remainingFeatures,
                                    setprefix = 'Remaining',
                                    featureColNames = featureColNames)
                                    
##
featureColNames <- colnames(virusFeatures)[grepl(ALL_FEATURES, colnames(virusFeatures))]
featureColNames <- featureColNames[! featureColNames %in% c('N', 'ATG.Bias')] # These don't vary

virusFeatures <- virusFeatures %>%
    select(UniversalName, Strain, featureColNames)

feature_f <- paste(calculated_dir, "GenomicFeatures-Virus.rds", sep="/")
saveRDS(virusFeatures, feature_f)

isgFeatures$GeneSet <- 'ISG'
housekeepingFeatures$GeneSet <- 'Housekeeping'
remainingFeatures$GeneSet <- 'Remaining'

humanFeatures <- bind_rows(isgFeatures, housekeepingFeatures, remainingFeatures) %>%
    select(-.data$Condition)

feature_f <- paste(calculated_dir, "GenomicFeatures-HumanCombined.rds", sep="/")
saveRDS(humanFeatures, feature_f)

featureDists <- isgDists %>%
    full_join(housekeepingDists, by = c('UniversalName', 'Strain')) %>%
    full_join(remainingDists, by = c('UniversalName', 'Strain'))

feature_f <- paste(calculated_dir, "GenomicFeatures-Distances.rds", sep="/")
saveRDS(featureDists, feature_f)