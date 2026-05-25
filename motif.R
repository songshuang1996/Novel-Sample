# =============================================================================
# motif.R — Motif filtering, merging, and logo plot generation
# =============================================================================
# PURPOSE:
#   This script is called automatically by the pipeline (via demo.ipynb) after
#   MEME Suite has finished running. It can operate in two modes:
#
#   MODE 1 — Single-folder mode:
#     If a file called "meme.txt" exists in the current directory, the script
#     reads it directly, generates a logo plot, saves it as "novel_motif.png",
#     and exits. This is used when set_length == 0 (comparative-only mode).
#
#   MODE 2 — Multi-folder mode:
#     If no "meme.txt" exists in the current directory, the script iterates
#     over every sub-folder (each representing one gene combination), reads
#     its MEME result, combines all motifs, filters for quality, merges
#     similar ones, and produces a final summarised logo plot.
#
# INPUT:
#   Current working directory must contain either:
#     - meme.txt                  (MODE 1)
#     - <gene_combo>/meme.txt     (MODE 2, one sub-folder per combination)
#
# OUTPUT:
#   novel_motif.png   — logo plot of all significant (filtered/merged) motifs
#   Motif_<n>         — frequency tables for each motif (MODE 2 only)
#   summary.txt       — combined frequency summary (written by Python wrapper)
# =============================================================================

library(universalmotif)   # motif I/O, filtering, merging, and visualization
library(ggplot2)          # backend for ggsave()
library(stringr)          # string manipulation for motif name parsing

# ── Clean workspace ────────────────────────────────────────────────────────────
rm(list = ls())

# ── Collect items in the current directory ─────────────────────────────────────
items <- list.files()

# =============================================================================
# MODE 1: Single meme.txt in the current directory (comparative-only run)
# =============================================================================
if ("meme.txt" %in% items) {
  message("MODE 1: Found meme.txt in current directory — single-run mode.")

  # Read all motifs from the MEME output file
  motif <- read_meme("meme.txt")

  # Generate and save the logo plot
  # Height scales with the number of motifs so each logo has enough space
  logplot <- view_motifs(motif)
  ggsave(
    "novel_motif.png",
    plot      = logplot,
    width     = 50,
    height    = length(motif) * 8,
    units     = "cm",
    device    = "png",
    limitsize = FALSE
  )

  message("Saved novel_motif.png — done.")
  stop()   # exit; nothing more to do in single-run mode
}

# =============================================================================
# MODE 2: Multiple sub-folders, each with its own MEME result
# =============================================================================
message("MODE 2: Scanning sub-folders for MEME results...")

motiftogether <- c()   # accumulator for all motifs across all combinations

for (folder in items) {

  # Each item should be a directory containing a meme.txt
  meme_path <- sprintf("%s/meme.txt", folder)

  # Read the MEME output for this gene combination
  motif <- tryCatch(
    read_meme(meme_path),
    error = function(e) {
      message(sprintf("  Skipping '%s': could not read meme.txt (%s)", folder, e$message))
      return(list())
    }
  )

  # Skip folders with no significant motifs
  if (length(motif) == 0) next

  # Ensure motif is always a list (universalmotif returns a single object when
  # there is only one motif, not a list)
  if (length(motif) == 1) motif <- c(motif)

  # Rename each motif to encode its origin folder and index.
  # Format: "<folder>-<total_motifs_in_folder>-<motif_index>"
  # This name is later parsed to reconstruct how many combinations contributed.
  for (i in seq_along(motif)) {
    motif[[i]]["name"] <- sprintf("%s-%s-%s", folder, length(motif), i)
  }

  # Append to the global collection
  motiftogether <- c(motiftogether, motif)
  message(sprintf("  Loaded %d motif(s) from '%s'", length(motif), folder))
}

message(sprintf("Total motifs collected before filtering: %d", length(motiftogether)))

# =============================================================================
# Quality filtering
# =============================================================================
# Keep only motifs that are:
#   - At least 14 bp wide  (removes very short, low-specificity motifs)
#   - IC score >= 15 bits  (information content; ensures the motif is specific)
motiftogether <- filter_motifs(motiftogether, width = 14, icscore = 15)

# Trim low-information positions from the edges of each motif
# (positions with IC < 0.6 bits per position are removed)
motiftogether <- trim_motifs(motiftogether, min.ic = 0.6)

# Re-apply the width/IC filter after trimming, because trimming can shorten motifs
motiftogether <- filter_motifs(motiftogether, width = 14, icscore = 15)

message(sprintf("Motifs after quality filtering: %d", length(motiftogether)))

# =============================================================================
# Merge similar motifs
# =============================================================================
# Motifs with a similarity score >= 0.6 (Pearson correlation of PWM columns)
# are considered redundant and merged into a consensus motif.
# nthreads controls parallelism in the similarity computation.
motiftogether1 <- merge_similar(motiftogether, threshold = 0.6, nthreads = 64)

# Trim again after merging, as consensus motifs may have ragged edges
motiftogether1 <- trim_motifs(motiftogether1, min.ic = 0.6)

message(sprintf("Motifs after merging similar: %d", length(motiftogether1)))

# =============================================================================
# Annotate motifs and write per-motif frequency tables
# =============================================================================
# For each merged motif we reconstruct:
#   - How many distinct gene combinations contributed ("frequency of occurrence")
#   - Which proteins appeared most often across those combinations

count_motif  <- 0
motiftogether2 <- motiftogether1   # copy that will receive renamed motifs

process_motif <- function(motif_name, motif_index) {
  # Strip the "-<n_motifs>-<index>" suffix to recover the folder names
  protein_id_all <- str_replace_all(
    string  = motif_name,
    pattern = "-\\d+-\\d+",
    replacement = ""
  )

  # Count occurrences: merged names are joined by "/" (one per source combination)
  occur_times <- str_count(motif_name, "/") + 1

  # Split the compound name into individual protein IDs
  split_data <- str_split(protein_id_all, "/|-")[[1]]

  # Build a frequency table of how often each protein appears across combinations
  count_data <- sort(table(split_data), decreasing = TRUE)

  # Save the frequency table as a plain CSV
  write.table(
    count_data,
    file      = sprintf("Motif_%d", motif_index),
    row.names = FALSE,
    col.names = FALSE,
    sep       = ","
  )

  # Return the human-readable label for this motif
  sprintf("Motif %d — Frequency of occurrence: %d", motif_index, occur_times)
}

# Handle both single-motif and multi-motif cases (universalmotif API differs)
if (length(motiftogether1) > 1) {
  for (i in seq_along(motiftogether1)) {
    count_motif <- count_motif + 1
    motif_name  <- motiftogether1[[i]]["name"]
    label       <- process_motif(motif_name, count_motif)
    motiftogether2[[count_motif]]["name"] <- label
  }
} else {
  # Single merged motif — universalmotif returns a non-list object
  motif_name <- motiftogether1["name"]
  label      <- process_motif(motif_name, 1)
  motiftogether2["name"] <- label
}

# =============================================================================
# Generate the final logo plot
# =============================================================================
# text.size = 4 keeps motif names readable even when many motifs are stacked.
# Height scales with the number of motifs (120 px per motif).
logplot <- view_motifs(motiftogether2, text.size = 4)

ggsave(
  "novel_motif.png",
  plot      = logplot,
  width     = 1400,
  height    = length(motiftogether2) * 120,
  units     = "px",
  device    = "png",
  limitsize = FALSE
)

message(sprintf("Saved novel_motif.png with %d motif(s).", length(motiftogether2)))
