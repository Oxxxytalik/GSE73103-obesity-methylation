# =====================================================================
# install_packages.R: install the package versions of the reference run
# (R 4.5.1, Windows 11 x64, Bioconductor 3.22; see session_info/ for the full sessionInfo()).
#
# Usage: Rscript install_packages.R
# =====================================================================

if (getRversion() < "4.5.0")
  warning("The reference run used R 4.5.1; Bioconductor 3.22 requires R 4.5.x.", call. = FALSE)
repos <- "https://cloud.r-project.org"
for (p in c("remotes", "BiocManager")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}

# ---- CRAN packages (exact versions) ----
cran <- c(
  data.table  = "1.18.2.1",
  ggplot2     = "4.0.2",
  patchwork   = "1.3.2",
  matrixStats = "1.5.0",
  uwot        = "0.2.5",
  lmtest      = "0.9-40",
  ggrepel     = "0.9.8",
  xml2        = "1.5.2",
  httr        = "1.4.8",
  glmnet      = "5.1",
  BiasedUrn   = "2.0.12",   # needed by missMethyl::gometh
  gPCA        = "1.0"       # archived on CRAN, installed from the archive
)
for (p in names(cran)) {
  have <- if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA
  if (!identical(have, cran[[p]])) {
    remotes::install_version(p, version = cran[[p]], repos = repos, upgrade = "never")
  }
}

# ---- R Markdown (any recent version; they only render the reports and do not change results) ----
for (p in c("knitr", "rmarkdown")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}

# ---- Bioconductor 3.22 ----
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
bioc <- c(
  EpiDISH = "2.26.0",
  minfi = "1.56.0",
  IlluminaHumanMethylation450kanno.ilmn12.hg19 = "0.6.1",
  qvalue = "2.42.0",
  missMethyl = "1.44.0",
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19 = "0.6.0",   # loaded by missMethyl
  org.Hs.eg.db = "3.22.0",
  GO.db = "3.22.0"
)
BiocManager::install(names(bioc), ask = FALSE, update = FALSE)

# ---- Report: installed vs expected versions ----
want <- c(cran, bioc)
got <- character(length(want))
for (i in seq_along(want)) {
  p <- names(want)[i]
  got[i] <- if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA
}
check <- data.frame(package = names(want), expected = unname(want), installed = got,
                    match = got == unname(want))
print(check, row.names = FALSE)
if (!all(check$match %in% TRUE)) message("Some versions differ from the reference run (see the table).")
