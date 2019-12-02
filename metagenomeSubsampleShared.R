#!/usr/bin/env Rscript
# metagenomeSubsampleShared.R
# William L. Close
# Pat Schloss Lab
# University of Michigan


# Setting environment -----------------------------------------------------

# Parsing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Variables defined by user
sharedFile <- args[1]

if(!file.exists(sharedFile)){
	stop(paste(sharedFile, "does not exist"))
}

nReads <- suppressWarnings(as.integer(args[2])) # Number of reads to subsample to

# if(is.na(nReads)){
# 	stop("The number of reads was not provided")
# }
# or...
if(is.na(nReads)){
	print("A numeric value was not given for the nReads argument so we will select the size of the smallest sample")
}

# Other variables
outDir <- "./"



# Loading dependencies ----------------------------------------------------

if(!require(tidyverse)){
	install.packges("tidyverse")
	library(tidyverse)
}



# Functions ---------------------------------------------------------------

# Function for removing and reporting samples without enough reads before subsampling
filter_shared_by_reads <- function(shared_df, n_reads_int) {

  # Finding smallest number of sample reads above specified subsampling threshold
  subsample_reads <- shared_df %>%
    filter(numReads >= n_reads_int) %>%
    pull(numReads) %>%
    min()

  # Pulling sample names WITHOUT enough reads for subsampling
  remove_samples <- shared_df %>%
    filter(numReads < subsample_reads) %>% # Finding samples that don't have enough reads for subsample
    pull(sample)

  # Removing samples without enough reads from downstream analysis
  filtered_shared <- shared_df %>%
    filter(!(sample %in% remove_samples))

	if(length(remove_samples) > 0){
		print(paste("Warning: these samples were removed because they had fewer than", n_reads_int, "reads:", paste(remove_samples, collapse=" ")))
	}

  return(filtered_shared)

}


# Function for subsampling a specific row (sample) of shared file
subsample_shared_row <- function(shared_df, sample_chr) {

	print(paste0("   Processing sample ", sample_chr, "..."))

  # Pulling out list of bins and sample counts for creating sampling vector
  bin_names <- str_subset(names(shared_df), "^Bin\\d+$") # Pulling bin names from col names
  sample_counts <- as.numeric(shared_df[shared_df$sample == sample_chr, bin_names]) # Pulling bin counts for one row at a time
  subsample_reads <- min(shared_df$numReads)

  # Creating sampling vector
  sampling_vector <- rep(bin_names, sample_counts) # Repeats each bin name for as many read counts for that bin

  # Subsampling sampling vector to n_reads
  subsampled_row <- sample(x = sampling_vector, size = subsample_reads) %>% # Subsampling to desired read depth
    table() %>% # Counting instances of each bin
    enframe(value = sample_chr)

  return(subsampled_row)

}


# Function for subsampling entire shared file to desired read depth
subsample_shared <- function(shared_df) {

  # Subsampling each row and joining into final subsampled shared file
  subsampled_shared <- map(as.character(shared_df$sample), subsample_shared_row, shared_df = shared_df) %>%
    reduce(full_join, by = "name") %>% # Combines all of the subsampled sample dfs keepin all bins with counts from all samples
    right_join(tibble(name = str_subset(names(shared_df), "^Bin\\d+$")), by = "name") %>% # Joining to list of all bins to fill in missing bins
    gather(key = "sample", value = "count", -name) %>% # Transforming to long df
    mutate(count = replace_na(count, 0), # Changing NAs from join to 0
           count = as.integer(count)) %>% # Recoding counts as integers
    mutate(numBins = length(unique(name))) %>% # Adding column of total bin counts
    spread(key = name, value = count) # Spreading data back out to shared format

  return(subsampled_shared)

}



# Analysis ----------------------------------------------------------------

# Creating output directory if it doesn't already exist
dir.create(outDir, recursive = TRUE, showWarnings=FALSE) # more platform independent

# Reading in data files
print("Reading in shared file...")
shared <- read_tsv(sharedFile, col_types = cols())


# This will select the size of the smallest library if the size was not provided
if(is.na(nReads)){
	nReads <- min(shared$numReads)
	print(paste("Samples subsampled to", nReads, " reads per sample"))
}


# Removing samples without enough reads
print("Removing samples without enough reads...")
filtered_shared <- filter_shared_by_reads(shared, nReads)

# Subsampling the shared file
print("Subsampling the shared file...")
sub_shared <- subsample_shared(filtered_shared)

# Write out the subsampled shared file
print("Write out the subsampled shared file")
write_tsv(sub_shared, path = paste0(outDir, "metagenome.subsample.shared"))
