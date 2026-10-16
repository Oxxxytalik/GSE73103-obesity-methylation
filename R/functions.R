# =====================================================================
# R/functions.R: packages, file names, colours, plot style and helper functions
#
# Every report in analysis/ sources this file in its first code chunk.
# =====================================================================


# ---- Packages ----------------------------------------------------------
# Install them once with install_packages.R (exact versions of the reference run).
suppressPackageStartupMessages({
  library(data.table)    # fread(): fast reading of the 2.9 GB beta-value file
  library(matrixStats)   # rowVars(), colMaxs()
  library(ggplot2)       # all figures
  library(patchwork)     # several ggplots in one figure (p1 | p2, p1 / p2)
  library(ggrepel)       # text labels that do not overlap
  library(lmtest)        # lrtest(): likelihood-ratio test of two nested models
  library(qvalue)        # Storey q-values
  library(uwot)          # UMAP
  library(gPCA)          # guided PCA test
  library(EpiDISH)       # blood cell composition from methylation
  library(glmnet)        # elastic-net logistic regression
})


# ---- Folders and input files -------------------------------------------
# All paths are relative to the repository root (the reports set it as the working directory).
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
result_path <- function(file_name) file.path("results", file_name)

# Folder with the three input files. By default the repository root; if the beta file is not there,
# the folders one and two levels above are tried. Any other folder can be set before running:
#   Sys.setenv(GSE73103_DATA_DIR = "C:/path/to/data")
DATA_DIR <- Sys.getenv("GSE73103_DATA_DIR", unset = "")
if (DATA_DIR == "") {
  DATA_DIR <- "."
  for (candidate in c(".", "..", "../..")) {
    if (file.exists(file.path(candidate, "GSE73103_beta.txt"))) {
      DATA_DIR <- candidate
      break
    }
  }
}
BETA_FILE  <- file.path(DATA_DIR, "GSE73103_beta.txt")                # GEO supplementary file: beta values of 355 samples
PHENO_FILE <- file.path(DATA_DIR, "GSE73103_phenotype_genotype.csv")  # phenotypes + 52 SNP genotypes (downloaded if missing)
SNP_DOCX   <- file.path(DATA_DIR, "13073_2015_225_MOESM1_ESM.docx")   # Additional file 1 of Voisin et al. 2015 (the 52 SNPs)
SNP_TABLE  <- result_path("snp_table.csv")        # made from SNP_DOCX by read_snp_table_docx()
BETA_CACHE <- result_path("beta_matrix.rds")      # the beta matrix saved in R format (fast to reload)


# ---- Colours (Okabe-Ito palette, safe for colour-blind readers) ---------
PAL_WEIGHT <- c("Normal-weight" = "#0072B2", "Overweight" = "#E69F00", "Obese" = "#D55E00")
PAL_SEX    <- c("#CC79A7", "#009E73")                 # Female, Male (order of the Sex factor levels)
PAL_FLAG   <- c("FALSE" = "grey55", "TRUE" = "#D55E00")
PAL_GENO   <- c("0" = "#9ECAE1", "1" = "#4292C6", "2" = "#08306B")   # 0/1/2 copies of an allele


# ---- Plot style -----------------------------------------------------------
theme_gse <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "right")
}

# The same style with the title aligned to the left edge of the whole figure.
theme_left <- function() {
  theme_gse() + theme(plot.title.position = "plot")
}

# Break a long text into lines of at most `width` characters (used for subtitles).
wrap_text <- function(x, width) {
  paste(strwrap(x, width), collapse = "\n")
}

# Save a figure as PNG (300 dpi) and PDF in figures/.
# `p` is a ggplot (or patchwork) object, or a function that draws a base-R plot.
save_fig <- function(p, name, width = 6.5, height = 5, dpi = 300) {
  file_png <- file.path("figures", paste0(name, ".png"))
  file_pdf <- file.path("figures", paste0(name, ".pdf"))
  if (inherits(p, "ggplot")) {
    # a long subtitle of a single plot is wrapped to the figure width (otherwise it is cut off)
    if (!inherits(p, "patchwork") && is.character(p$labels$subtitle)) {
      p$labels$subtitle <- wrap_text(p$labels$subtitle, floor(width * 11))
    }
    ggsave(file_png, p, width = width, height = height, dpi = dpi, bg = "white")
    ggsave(file_pdf, p, width = width, height = height)
  } else if (is.function(p)) {
    png(file_png, width = width, height = height, units = "in", res = dpi); p(); dev.off()
    pdf(file_pdf, width = width, height = height); p(); dev.off()
  } else {
    stop("save_fig: expected a ggplot object or a function that draws a plot")
  }
  cat("Figure saved:", file_png, "\n")
  invisible(file_png)
}

# Legend labels with the group size, e.g. "Normal-weight (n=268)".
level_labels_n <- function(f) {
  counts <- table(f)
  setNames(sprintf("%s (n=%d)", names(counts), as.integer(counts)), names(counts))
}

scale_colour_weight <- function(f) {
  labels <- level_labels_n(f)
  scale_colour_manual(values = PAL_WEIGHT[names(labels)], breaks = names(labels),
                      labels = unname(labels), name = "Weight category", drop = FALSE)
}

scale_colour_sex <- function(f) {
  labels <- level_labels_n(f)
  colours <- setNames(PAL_SEX[seq_along(labels)], names(labels))
  scale_colour_manual(values = colours, breaks = names(labels), labels = unname(labels),
                      name = "Sex", drop = FALSE)
}

# Axis label of a principal component with its share of the variance.
pc_label <- function(pca, k) {
  percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  sprintf("PC%d (%.1f%% of variance)", k, percent[k])
}

# Sort rows so that the rarest group is drawn last (on top):
# otherwise the 17 Obese points are hidden under the 268 Normal-weight points.
order_rare_last <- function(df, f) {
  group_size <- as.integer(table(f)[as.character(f)])
  df[order(-group_size), , drop = FALSE]
}


# ---- Small statistics helpers --------------------------------------------

# Genomic inflation factor lambda: median chi-square statistic of the observed p-values
# divided by the median expected under the null, qchisq(0.5, 1) = 0.455. lambda = 1 means no inflation.
lambda_gc <- function(p) {
  p <- p[!is.na(p) & p > 0]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

# Area under the ROC curve from ranks (Mann-Whitney formula).
# y = 0/1 labels, score = higher means "more likely 1". 0.5 = chance, 1 = perfect.
auc_value <- function(score, y) {
  r <- rank(score)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# First (up to 3) distinct gene names of the 450K annotation column UCSC_RefGene_Name.
short_gene_names <- function(refgene) {
  parts <- strsplit(as.character(refgene), ";")
  out <- character(length(parts))
  for (i in seq_along(parts)) {
    genes <- parts[[i]]
    genes <- unique(genes[!is.na(genes) & genes != ""])
    if (length(genes) > 0) out[i] <- paste(head(genes, 3), collapse = ";")
  }
  out
}


# ---- Fast linear regression over all CpGs ---------------------------------
# For every CpG (row of M) this gives the same numbers as
#     summary(lm(M[cpg, ] ~ y + covariates))   for the coefficient of y
# (estimate, standard error, t and p), but for 397,615 CpGs in seconds instead of hours.
# How: the covariates and the intercept are first removed from y and from every CpG
# (projection on the QR basis Q of the covariate matrix). After that, the slope of y is
#     b = sum(y_res * M_res) / sum(y_res^2)              (Frisch-Waugh-Lovell theorem)
# and the residual variance is (RSS of the model without y - b^2 * sum(y_res^2)) / df.
# Each report checks the result against lm() on a few random CpGs.
#
# Arguments
#   M          CpG x sample matrix (M-values)
#   y          variable of interest (e.g. 0 = Normal-weight, 1 = Overweight/Obese)
#   covariates data frame of adjustment variables (one row per sample)
#   keep       TRUE/FALSE per sample: analyse a subset (e.g. women only)
#   y_perm     optional matrix of permuted versions of y (one column per permutation);
#              they are tested in the same pass, for permutation thresholds
#   keep_perm_t TRUE = return the full CpG x permutation matrix of t-statistics
# Covariates that are constant in the subset (e.g. Sex among women only) are dropped.
ols_scan <- function(M, y, covariates, keep = rep(TRUE, ncol(M)), y_perm = NULL,
                     keep_perm_t = FALSE, chunk_size = 40000) {
  X <- cbind(1, as.matrix(covariates[keep, , drop = FALSE]))
  sd_cols <- apply(X[, -1, drop = FALSE], 2, sd)
  X <- X[, c(TRUE, sd_cols > 0), drop = FALSE]
  Q <- qr.Q(qr(X))
  df <- sum(keep) - ncol(X) - 1                     # residual degrees of freedom

  # column 1 = real y, further columns = permuted y; covariates removed from all of them
  Y_all <- cbind(y[keep], y_perm[keep, , drop = FALSE])
  G <- Y_all - Q %*% crossprod(Q, Y_all)
  ss_g <- colSums(G^2)
  n_perm <- ncol(G) - 1

  n_cpg <- nrow(M)
  est <- rep(NA_real_, n_cpg)
  se <- rep(NA_real_, n_cpg)
  tstat <- rep(NA_real_, n_cpg)
  max_abs_t_perm <- rep(0, n_perm)                  # largest |t| over all CpGs, per permutation
  t_perm <- NULL
  if (keep_perm_t) t_perm <- matrix(NA_real_, n_cpg, n_perm)

  # the CpGs are processed in chunks of 40,000 rows to limit memory use
  for (first in seq(1, n_cpg, by = chunk_size)) {
    rows <- first:min(first + chunk_size - 1, n_cpg)
    Y <- M[rows, keep, drop = FALSE]
    Y_res <- Y - (Y %*% Q) %*% t(Q)                 # methylation with covariates removed
    rss0 <- rowSums(Y_res^2)
    b <- sweep(Y_res %*% G, 2, ss_g, "/")           # slope: divide column j by ss_g[j]
    rss <- rss0 - sweep(b^2, 2, ss_g, "*")
    s2 <- pmax(rss, 1e-12) / df                     # residual variance
    tt <- b / sqrt(sweep(s2, 2, ss_g, "/"))         # t = b / SE, with SE = sqrt(s2 / ss_g)
    est[rows] <- b[, 1]
    se[rows] <- sqrt(s2[, 1] / ss_g[1])
    tstat[rows] <- tt[, 1]
    if (n_perm > 0) {
      max_abs_t_perm <- pmax(max_abs_t_perm, colMaxs(abs(tt[, -1, drop = FALSE])))
      if (keep_perm_t) t_perm[rows, ] <- tt[, -1]
    }
  }
  list(est = est, se = se, t = tstat, p = 2 * pt(-abs(tstat), df), df = df,
       max_abs_t_perm = max_abs_t_perm, t_perm = t_perm)
}

# Compare ols_scan() with lm() on a few CpGs (rows). Stops if they disagree.
check_scan_with_lm <- function(scan, M, y, covariates, rows) {
  d <- cbind(covariates, y_of_interest = y)
  diff_beta <- numeric(0)
  diff_log10p <- numeric(0)
  for (i in rows) {
    fit <- summary(lm(as.numeric(M[i, ]) ~ ., data = d))$coefficients["y_of_interest", ]
    diff_beta <- c(diff_beta, abs(fit[["Estimate"]] - scan$est[i]))
    diff_log10p <- c(diff_log10p, abs(log10(fit[["Pr(>|t|)"]]) - log10(scan$p[i])))
  }
  cat(sprintf("Check against lm() on %d CpGs: max |difference in estimate| = %.1e, max |difference in log10 p| = %.1e\n",
              length(rows), max(diff_beta), max(diff_log10p)))
  stopifnot(max(diff_beta) < 1e-6, max(diff_log10p) < 1e-4)
}


# ---- Check that a previous report has run in this R session ---------------
need_objects <- function(names, report) {
  missing <- names[!vapply(names, exists, logical(1))]
  if (length(missing) > 0) {
    stop("Objects from ", report, " are not in memory (", paste(missing, collapse = ", "), "). ",
         "Render the reports in order in one R session: source(\"run_all.R\").", call. = FALSE)
  }
}
