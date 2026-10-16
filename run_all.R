# =====================================================================
# run_all.R: render all reports in order, in one R session
#
# Each report (analysis/*.Rmd) is rendered to a Markdown file next to it (analysis/*.md)
# The code shows it's output and the figures. Later reports use objects created by
# earlier ones (pheno, mvals_full, dosage_mat, ...), so they are rendered one after another
# in the same session: envir = globalenv() keeps every object in memory for the next report.
#
# Run time of the reference run: about 1 hour (most of it in report 6, the 100 permutations
# of the prediction pipeline).
# =====================================================================

# ---- Paths on the computer of the reference run (setwd() for correct results on the other devise) ----
for (folder in c(".",
                 "C:/Users/User/Downloads/GSE73103-obesity-methylation",
                 "C:/Users/User/Downloads/GSE73103-obesity-methylation/GSE73103-obesity-methylation")) {
  if (file.exists(file.path(folder, "analysis/01_data_pca_umap.Rmd"))) {
    setwd(folder)
    break
  }
}
# Folder with GSE73103_beta.txt, GSE73103_phenotype_genotype.csv and 13073_2015_225_MOESM1_ESM.docx
if (file.exists("C:/Users/User/Downloads/GSE73103_beta.txt")) Sys.setenv(GSE73103_DATA_DIR = "C:/Users/User/Downloads")

if (!file.exists("analysis/01_data_pca_umap.Rmd"))
  stop("Run run_all.R from the repository root (the folder that contains analysis/).", call. = FALSE)
cat("Repository:", getwd(), "\n")

reports <- c(
  "analysis/01_data_pca_umap.Rmd",     # data, QC, PCA, outliers, gPCA, cell composition, UMAP
  "analysis/02_meqtl.Rmd",             # cis-meQTL, genotype permutations, GoDMC positive control
  "analysis/03_coherence.Rmd",         # genotype-free coherence check of the paper's SNP-CpG pairs
  "analysis/04_ewas.Rmd",              # EWAS Overweight/Obese vs Normal-weight
  "analysis/05_ewas_followups.Rmd",    # known BMI CpGs, obesity loci, DMR, sensitivity, GO
  "analysis/06_power_prediction.Rmd"   # power, covariate-adjusted prediction, summary table
)

for (report in reports) {
  start <- Sys.time()
  message("\n######## ", report, "  (started ", format(start, "%H:%M:%S"), ")")
  rmarkdown::render(report, envir = globalenv(), quiet = TRUE)
  message(sprintf("######## %s finished in %.1f min", report,
                  as.numeric(difftime(Sys.time(), start, units = "mins"))))
}

writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("\nDone. Reports: analysis/*.md, figures: figures/, tables: results/.")
