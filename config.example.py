# =============================================================================
# config.example.py  —  Path configuration template
# =============================================================================
# INSTRUCTIONS:
#   1. Copy this file to config.py:
#          cp config.example.py config.py
#   2. Fill in all values below to match your system.
#   3. config.py is git-ignored and will not be committed.
# =============================================================================

# Root directory where all pipeline data will be stored.
# The pipeline will create these sub-folders inside BASE_DIR automatically:
#   0_gbk/          — downloaded GenBank genome files
#   0_protein/      — extracted protein FASTA files
#   0_promoter/     — extracted promoter FASTA files
#   99_User_project/ — per-job working directories
BASE_DIR = "/your/working/directory"

# Full path to the MEME executable (not just the bin folder).
# Example: "/opt/conda/envs/novel-cis/bin/meme"
MEME_BIN = "/path/to/bin/meme"

# Full path to the Rscript executable.
# Example: "/opt/conda/envs/R/bin/Rscript"
RSCRIPT = "/path/to/Rscript"

# Full path to the motif.R script in this repository.
# Example: "/home/user/Novel-Sample/motif.R"
MOTIF_R = "/path/to/Novel-Sample/motif.R"

# Email address used for NCBI Entrez API calls (required by NCBI).
EMAIL = "your@email.com"

# Number of CPU threads for parallel steps (DIAMOND, CD-HIT, MEME, multiprocessing).
# Set to the number of cores available on your machine.
THREADS = 32

# Port for the PyWebIO web server.
PORT = 8848
