# =============================================================================
# config.R — GWASphere configuration
# =============================================================================
# Set GWASPHERE_BASE_DIR to point to the app folder.
# Default: $HOME/gwasphere
# =============================================================================

BASE <- Sys.getenv("GWASPHERE_BASE_DIR", file.path(Sys.getenv("HOME"), "gwasphere"))

REGENIE_APP_DIR <- file.path(BASE, "1.Regenie")
VIEWER_APP_DIR  <- file.path(BASE, "3.Visualization_webapp")
PIPELINE_SCRIPT <- file.path(BASE, "2.Downstream", "run_pipeline.sh")
