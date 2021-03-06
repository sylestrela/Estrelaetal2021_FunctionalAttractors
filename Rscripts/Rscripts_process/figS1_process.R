" R script to do some preprocessing for figure S1 in 
Estrela et al (2021) Functional attractors in microbial community assembly."

## ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(reshape2)
library(tidyr)
## ---------------------------------------------------------------------------

#.. Assign path to data
#. usern = path to user directory

path.d<-paste0(usern, "/data/data_main")
#.. Assign path to plots output
path.p<- paste0(usern,"/Rplots")
#.. Assign path to data figures
path.ds<- paste0(usern,"/data/data_figures")
#.. Assign path to data processed
path.dp<- paste0(usern,"/data/data_processed")

#.. Read data 
setwd(path.d)

# Match isolates' 16S sequences to community ESV
isolates_ID_RDP <- fread("Estrela_2021_isolates_16S_taxa.csv")
communities_abundance_syl <- fread("Goldford_2018_16S_data_table_glucose_RelAbund_ALL.csv")

# Modify column name to match function needs
isolates_ID_RDP <- isolates_ID_RDP %>%
    mutate(Community = InocRep, ID = SangerID)
communities_abundance_syl <- communities_abundance_syl %>%
    mutate(Community = InocRep, RelativeAbundance = Relative_Abundance, CommunityESVID = ESV_ID)

# Match isolate Sanger seq to community ESV ----
## R function for pairwise alignment and counting bp gap and mismatch
alignment_bp_count <- function(seq1, seq2, type = "local"){
    # Make pairwise alignment
    algn <- Biostrings::pairwiseAlignment(seq1, seq2, type = type)

    # Make a df
    temp_df <- tibble(
        s1 = strsplit(as.character(Biostrings::pattern(algn)), '')[[1]],
        s2 = strsplit(as.character(Biostrings::subject(algn)), '')[[1]]
    )

    # Find difference; ignore N (which means any nucleotide)
    temp_df <- temp_df %>% mutate(Diff = ((s1 != s2) & s1 != 'N' & s2 != 'N'))

    # Find gaps
    ngap <- sum(temp_df$s1 == "-" | temp_df$s2 == "-")

    # Find mismatch
    nmis <- sum((temp_df$s1 != temp_df$s2) & temp_df$s1 != "-" & temp_df$s2 != "-" & temp_df$s1 != "N" & temp_df$s2 != "N")

    # Length of consensus
    nlength <- nrow(temp_df)

    #
    c(ConsensusLength = nlength, BasePairGap = ngap, BasePairMismatch = nmis, AlignmentScore = algn@score) %>%
        return()
}
## Example code for match one pair of sequences
# seq1 <- isolates_ID_RDP$Sequence[1]
# seq2 <- isolates_ID_RDP$Sequence[3]
# alignment_bp_count(seq1, seq2) # R function from `sequence_pairwise_alignment.R`
# printPairwiseAlignment(Biostrings::pairwiseAlignment(seq1, seq2))

## R function for align isolates' Sanger sequences to ESVs within communities, based on different alignment methods
align_sequence <- function(communities_abundance_align, isolates_ID_RDP, alignment_type = "global") {
    # Community names. All communties with isolates
    comm_name <- unique(isolates_ID_RDP$Community)

    # Empty list for alignment
    sequences_alignment_list_func <- rep(list(NA), length(comm_name)) %>% setNames(comm_name)

    # Loop through comunities
    for (i in 1:length(comm_name)) {
        # Subset for community ESV
        temp_comm <- communities_abundance_align %>%
            filter(Community == comm_name[i]) %>%
            select(Community, Family, RelativeAbundance, CommunityESVID, ESV)

        # Isolate sequences
        temp_iso <- isolates_ID_RDP %>%
            filter(Community == comm_name[i]) %>%
            # Change variable name used in isolates to avoid duplicate variable name
            mutate(IsolateID = ID, IsolateGenus = Genus, IsolateSangerSequence = Sequence) %>%
            select(Community, IsolateID, IsolateGenus, IsolateSangerSequence)

        # Merge two dfs: it will auto create df with nrow(temp_df) = nrow(temp_comm) * nrow(temp_iso)
        temp_df <- temp_comm %>% left_join(temp_iso, by = "Community")
        cat("\n", comm_name[i], "\n number of ESVs = ", nrow(temp_comm), "\n number of isolates = ", nrow(temp_iso), "\n pairs = ", nrow(temp_df), "\n")

        # Add five new columns to temp_df
        temp_df$AlignmentType <- "NA"
        temp_df[,c("ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <- NaN

        # Pairwise alignment
        for (j in 1:nrow(temp_df)) {
            # Skip isolates that don't have Sanger sequences
            if (temp_df$IsolateSangerSequence[j] == "" | is.na(temp_df$IsolateSangerSequence[j])) next

            # Alignment
            result <- alignment_bp_count(temp_df$ESV[j], temp_df$IsolateSangerSequence[j], type = alignment_type)
            temp_df[j,c("AlignmentType","ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <-
                tibble(AlignmentType = alignment_type,
                       ConsensusLength = result[1],
                       BasePairGap = result[2],
                       BasePairMismatch = result[3],
                       AlignmentScore = result[4])

            #
            if ((j %% 10) == 0) cat(j, " ")
        }

        # Assignment
        sequences_alignment_list_func[[i]] <- temp_df
        cat("\n")
    }

    # Merge list
    sequences_alignment <- rbindlist(sequences_alignment_list_func[!is.na(sequences_alignment_list_func)])  %>% as_tibble

    return(sequences_alignment)
}

## Align Sanger to ESVs
alignment_type <- c("global", "local", "overlap", "global-local", "local-global")
sequences_alignment_list <- rep(list(NA), length(alignment_type))
for (i in 1:length(alignment_type)) {
    sequences_alignment_list[[i]] <- align_sequence(communities_abundance_syl, isolates_ID_RDP, alignment_type = alignment_type[i])
    print(alignment_type[i])
}
## Merge list
sequences_alignment_syl <- rbindlist(sequences_alignment_list) %>% as_tibble
fwrite(sequences_alignment_syl, paste0(path.dp,"/figS1_sequences_alignment.csv"))

#' Filter for matched ESV-isolate and plot the relative abundances
communities_name_pool <- c(paste0("C", 1:12, "Rpool"), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))
families_name <- c("Aeromonadaceae", "Alcaligenaceae", "Bradyrhizobiaceae", "Brucellaceae", "Burkholderiaceae", "Caulobacteraceae", "Cellulomonadaceae", "Chitinophagaceae", "Chthoniobacteraceae", "Comamonadaceae", "Cryomorphaceae", "Enterobacteriaceae", "Enterococcaceae", "Flavobacteriaceae", "Hyphomicrobiaceae", "Listeriaceae", "Microbacteriaceae", "Moraxellaceae", "Nocardiaceae", "Obscuribacterales.17", "Oxalobacteraceae", "Paenibacillaceae", "Phyllobacteriaceae", "Porphyromonadaceae", "Pseudomonadaceae", "Rhizobiaceae", "Sanguibacteraceae", "Sphingobacteriaceae", "Sphingomonadaceae", "Xanthomonadaceae")

# Read data
## Sanger to ESV alignment
sequences_alignment_syl <- fread(paste0(path.dp,"/figS1_sequences_alignment.csv")) %>%
    mutate(Community = ordered(Community, levels = communities_name_pool))

# Find the match with highest alignment score
sequences_abundance_list <- rep(list(NA), 3)
allow_mismatch <- c(0:2, Inf)

for (i in 1:4) {
    sequences_abundance_list[[i]] <-
        sequences_alignment_syl %>%
        filter(AlignmentType == "local") %>%
        # Filter for BasePairMatch
        filter(BasePairMismatch <= allow_mismatch[i]) %>%
        # For each Sanger, find the Sanger-ESV match with highest alignment score
      dplyr::group_by(AlignmentType, IsolateID) %>%
      dplyr::arrange(desc(AlignmentScore)) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        # Remove duplicates that matches two Sangers to one ESV
      dplyr::group_by(AlignmentType, Community, CommunityESVID) %>%
        distinct(RelativeAbundance, .keep_all = T) %>%
        arrange(Community) %>%
        # Specify mismatch allowed
        mutate(AllowMismatch = allow_mismatch[i])
}

sequences_abundance <- rbindlist(sequences_abundance_list) %>% as_tibble

# Remove unessential variables
sequences_abundance <- sequences_abundance %>%
    mutate(AlignmentType = ordered(AlignmentType, levels = c("global", "local", "overlap", "global-local", "local-global"))) %>%
    select(AlignmentType, AllowMismatch, Community, IsolateGenus,
           RelativeAbundance, CommunityESVID, Family, ConsensusLength, BasePairGap, BasePairMismatch, AlignmentScore) %>%
    arrange(AlignmentType, AllowMismatch, Community)

fwrite(sequences_abundance, file = paste0(path.dp,"/figS1_sequences_abundance.csv"))


