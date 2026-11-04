

# ==========================================
# 1. Define & Install Required CRAN Packages
# ==========================================

# Added "remotes" to the list so GitHub installations work seamlessly
packages <- c("remotes", "nimble", "sf", "MCMCvis", "coda", "MASS", 
              "parallel", "dplyr", "ggplot2", "lubridate", "terra", "stars", "gstat")

# Function to install missing packages
install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Run the installation check
invisible(lapply(packages, install_if_missing))


# ==========================================
# 2. Install GitHub Packages (climateR & AOI)
# ==========================================
# These are not on CRAN and must be installed from GitHub.
#
# If installation fails with a rate-limit or authentication error, you may need a
# GitHub Personal Access Token. Create one at https://github.com/settings/tokens
# (no scopes required for public repositories) and register it with:
#
#     usethis::create_github_token()   # opens the page
#     gitcreds::gitcreds_set()         # stores the token
#
# See https://usethis.r-lib.org/articles/git-credentials.html for details.

remotes::install_github("mikejohnson51/AOI")
remotes::install_github("mikejohnson51/climateR")

# ==========================================
# 3. Load All Packages
# ==========================================

# Added suppressPackageStartupMessages to keep your console clean
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# Explicitly load the GitHub packages as well
suppressPackageStartupMessages(library(AOI))
suppressPackageStartupMessages(library(climateR))

