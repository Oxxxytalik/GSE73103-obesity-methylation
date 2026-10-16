# Blood DNA methylation, obesity-risk SNPs and weight: a re-analysis of GSE73103

The project analyses the public dataset **GEO [GSE73103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73103)** from Voisin et al. (2015, *Genome Medicine* 7:103):

- Illumina HumanMethylation450 whole-blood DNA methylation of **355 young people** (14–34 years), **397,615 CpGs** after GEO processing, no missing values;
- weight category of each person: 268 Normal-weight, 70 Overweight, 17 Obese (BMI itself is not available);
- genotypes of **52 obesity GWAS SNPs**.

**The code is written as **R Markdown reports** . 
---



## Repository layout

```
.
├── README.md
├── install_packages.R        # installs the package versions of the reference run
├── run_all.R                 # renders all reports in order, in one R session
├── R/
│   ├── functions.R           # packages, file names, colours, plot style, helper functions (incl. fast OLS)
│   └── input_files.R         # downloads the phenotype table; reads the SNP table from the paper's .docx
├── analysis/
│   ├── 01_data_pca_umap.Rmd      # data, QC, cohort, PCA, outliers, gPCA, cell composition, UMAP  -> figures 1–5
│   ├── 02_meqtl.Rmd              # cis-meQTL, genotype permutations, GoDMC positive control       -> figures 6, 7
│   ├── 03_coherence.Rmd          # genotype-free coherence check of the paper's SNP-CpG pairs     -> figure 8
│   ├── 04_ewas.Rmd               # EWAS Overweight/Obese vs Normal-weight                         -> figure 9
│   ├── 05_ewas_followups.Rmd     # known BMI CpGs, obesity loci, DMR scan, sensitivity, GO        -> figures 10, 11, 12, 13, 14
│   └── 06_power_prediction.Rmd   # power, covariate-adjusted prediction, summary table           -> figures 15, 16
├── data/
│   └── voisin2015_additional_file4_pairs.csv   # the paper's 107 significant SNP-CpG pairs (Additional file 4)
├── expected_results/         # result tables of the reference run, for comparison
└── session_info/             # reference run: input md5 checksums and run records (R version, RNG, sessionInfo)
```


## Data

| File | Source | md5 (reference run) | Size |
|---|---|---|---|
| `GSE73103_beta.txt` | GEO GSE73103 supplementary file with the processed β values: download from the series page, decompress, rename. 711 columns (probe ID, then `S1`, `Detection Pval S1`, … for 355 samples), 397,615 rows. | `b3ca0e92fdf9ea07659eca40869f0435` | 2,886,013,516 B |
| `GSE73103_phenotype_genotype.csv` | **Downloaded automatically** by report 1 if missing (`download_phenotypes()` reads the GEO page of each of the 355 samples). | `3718b22df92760ea6e61bbd35bc4acff` | 95,473 B |
| `13073_2015_225_MOESM1_ESM.docx` | Additional file 1 of Voisin et al. 2015 (the 52 SNPs with hg19 position and risk allele), from the article's supplementary material. | `dc31f6d47d21809524ad07c6f9fa181f` | 110,995 B |


Small reference tables are included in the code or in `data/`:

- `data/voisin2015_additional_file4_pairs.csv`: the 107 significant SNP–CpG pairs of Voisin et al. 2015, Additional file 4 (report 3);
- GoDMC blood meQTLs of rs6861681 and rs9939609, from <http://mqtldb.godmc.org.uk> (report 2);
- known BMI CpGs and their direction, from Wahl et al. 2017, Supplementary Table 3 (reports 4 and 5).

The 450K annotation (`IlluminaHumanMethylation450kanno.ilmn12.hg19`) and the EpiDISH blood reference (`centDHSbloodDMC.m`) come from their R packages.

---

## Software

Reference run: **R 4.5.1** (Windows 11 x64), **Bioconductor 3.22**, R's default random number generator (Mersenne-Twister / Inversion / Rejection).

| Package | Version | Used for |
|---|---|---|
| data.table | 1.18.2.1 | reading the β matrix (`fread`) |
| matrixStats | 1.5.0 | row variances |
| ggplot2 | 4.0.2 | all figures |
| patchwork | 1.3.2 | multi-panel figures |
| ggrepel | 0.9.8 | labels in the volcano plot |
| uwot | 0.2.5 | UMAP |
| gPCA | 1.0 | guided PCA test (archived on CRAN) |
| lmtest | 0.9-40 | likelihood-ratio tests (meQTL) |
| glmnet | 5.1 | elastic-net prediction |
| xml2 | 1.5.2 | reading the SNP table from the .docx |
| httr | 1.4.8 | downloading the GSM phenotype pages |
| BiasedUrn | 2.0.12 | needed by `missMethyl::gometh` |
| EpiDISH | 2.26.0 | blood cell-composition estimates |
| minfi | 1.56.0 | access to the 450K annotation |
| IlluminaHumanMethylation450kanno.ilmn12.hg19 | 0.6.1 | probe annotation (chromosome, position, gene) |
| qvalue | 2.42.0 | Storey q-values (meQTL) |
| missMethyl | 1.44.0 | GO enrichment corrected for probes per gene |
| org.Hs.eg.db / GO.db | 3.22.0 | GO annotation for missMethyl |
| IlluminaHumanMethylationEPICanno.ilm10b4.hg19 | 0.6.0 | loaded by missMethyl |
| knitr, rmarkdown | any recent | rendering the reports (no influence on the results) |

```r
Rscript install_packages.R   # installs these versions (CRAN archive + Bioconductor 3.22)
```

The full `sessionInfo()` of the reference run is in `session_info/run_record_*.txt`; a new run writes its own to `results/session_info.txt` and at the end of report 6.

---

## How to run
 Prior setwd() for the data files is recommended;

`run_all.R` renders the six reports **in order, in one R session** (`rmarkdown::render(..., envir = globalenv())`). Later reports use objects of earlier ones (`pheno`, `mvals_full`, `dosage_mat`, `ewas_summary`, …), so a single report cannot be knitted on its own in a fresh session; each report checks this and stops with a clear message.

- The first run reads the 2.9 GB text file and saves the matrix as `results/beta_matrix.rds` (about 1 GB); later runs load that file instead. Delete it to read the text file again.
- Run time of the reference run: about one hour, mostly the 100 permutations of the prediction pipeline (report 6). At least 16 GB of RAM is recommended.
- The main settings (thresholds, numbers of permutations, seeds) are at the top of each report in a `settings` chunk.

---

## Reproducibility

| Report | Seed | Used for |
|---|---|---|
| 01 | `SEED = 1234` | gPCA permutations, UMAP |
| 02 | `FU_SEED + 7 = 1241` | genotype permutations |
| 03 | `COH_SEED = 11` | random CpGs of the null band |
| 04 | `EW_SEED = 20261009` | label permutations, checks against `lm()`, volcano subsample |
| 05 | `GF_SEED = 20261010` | label permutations (DMR), check against `lm()`, random window sets |
| 06 | `GF2_SEED = 20261012` | check against `lm()`, ECDF subsample; the prediction uses `GF2_SEED + 100 + r`, `+ 200 + r`, `+ 300 + r` (repeat r) and `+ 1000 + i` (permutation i); the bootstrap CIs continue the random stream after `GF2_SEED + 303` |

---

## References

- Voisin S. et al. (2015) Many obesity-associated SNPs strongly associate with DNA methylation changes at proximal promoters and enhancers. *Genome Medicine* 7:103.
- Wahl S. et al. (2017) Epigenome-wide association study of body mass index, and the adverse outcomes of adiposity. *Nature* 541:81–86.
- Min J.L. et al. (2021) Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation (GoDMC). *Nature Genetics* 53:1311–1321.
- Reese S.E. et al. (2013) A new statistic for identifying batch effects in high-throughput genomic data that uses guided principal component analysis. *Bioinformatics* 29:2877–2883.
- Teschendorff A.E. et al. (2017) A comparison of reference-based algorithms for correcting cell-type heterogeneity in epigenome-wide association studies (EpiDISH). *BMC Bioinformatics* 18:105.
- Phipson B., Maksimovic J., Oshlack A. (2016) missMethyl: an R package for analyzing data from Illumina's HumanMethylation450 platform. *Bioinformatics* 32:286–288.
- Friedman J., Hastie T., Tibshirani R. (2010) Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software* 33:1–22.
- McInnes L., Healy J., Melville J. (2018) UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv:1802.03426.
- Cole T.J. et al. (2000) Establishing a standard definition for child overweight and obesity worldwide: international survey. *BMJ* 320:1240.
