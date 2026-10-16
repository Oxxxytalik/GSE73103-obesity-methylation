# =====================================================================
# R/input_files.R: create the two small input tables that are not plain downloads
#
#   download_phenotypes()  -> GSE73103_phenotype_genotype.csv
#       phenotypes (age, sex, weight category) and the genotypes of the 52 SNPs,
#       read from the GEO page of each of the 355 samples (GSM accessions)
#   read_snp_table_docx()  -> results/snp_table.csv
#       the paper's table of the 52 SNPs (position in hg19, risk allele, ...),
#       read from Additional file 1 of Voisin et al. 2015 (.docx)
#
# Both functions are called by analysis/01 and analysis/02 only if the output file is missing.
# =====================================================================


download_phenotypes <- function(out_file = PHENO_FILE) {
  # The 355 GSM accessions of the series, written as ranges (taken from the GEO series page)
  first <- c(1886364, 1886619, 1886713, 1886716, 1886720, 1886723, 1886726, 1886774)
  last  <- c(1886608, 1886693, 1886713, 1886716, 1886720, 1886723, 1886726, 1886803)
  gsm_ids <- character(0)
  for (k in seq_along(first)) gsm_ids <- c(gsm_ids, sprintf("GSM%d", first[k]:last[k]))
  stopifnot(length(gsm_ids) == 355, !anyDuplicated(gsm_ids))

  rows <- list()
  for (i in seq_along(gsm_ids)) {
    url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&view=brief&form=text",
                   gsm_ids[i])
    response <- httr::GET(url, httr::timeout(30))
    if (httr::status_code(response) != 200) stop(gsm_ids[i], ": HTTP ", httr::status_code(response))
    lines <- strsplit(httr::content(response, as = "text", encoding = "UTF-8"), "\n")[[1]]

    # lines look like "!Sample_characteristics_ch1 = key: value"
    char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
    key_value <- sub("^!Sample_characteristics_ch1 = ", "", char_lines)
    values <- trimws(sub("^[^:]+:", "", key_value))
    names(values) <- make.names(trimws(sub(":.*", "", key_value)))
    title <- sub("^!Sample_title = ", "", grep("^!Sample_title", lines, value = TRUE))

    rows[[i]] <- as.list(c(geo_accession = gsm_ids[i], title = title, values))
    if (i %% 25 == 0) cat(sprintf("...%d/%d samples\n", i, length(gsm_ids)))
    Sys.sleep(0.25)   # short pause between requests to the NCBI server
  }
  pheno_geno <- rbindlist(rows, fill = TRUE)
  stopifnot(nrow(pheno_geno) == 355)
  fwrite(pheno_geno, out_file)
  cat("Saved:", out_file, "\n")
}


read_snp_table_docx <- function(docx = SNP_DOCX, out_file = SNP_TABLE) {
  # A .docx file is a zip archive; the text of the document is in word/document.xml
  tmp <- tempfile("docx_")
  dir.create(tmp)
  utils::unzip(docx, files = "word/document.xml", exdir = tmp)
  doc <- xml2::read_xml(file.path(tmp, "word", "document.xml"))
  xml2::xml_ns_strip(doc)
  tables <- xml2::xml_find_all(doc, "//w:tbl")
  if (length(tables) != 1) stop("Expected exactly 1 table in the document, found ", length(tables))

  # text of one table cell: its paragraphs joined with line breaks
  cell_text <- function(cell) {
    paragraphs <- xml2::xml_find_all(cell, ".//w:p")
    texts <- character(0)
    for (p in paragraphs) texts <- c(texts, paste(xml2::xml_text(xml2::xml_find_all(p, ".//w:t")), collapse = ""))
    trimws(paste(texts, collapse = "\n"))
  }

  table_rows <- xml2::xml_find_all(tables[[1]], "./w:tr")
  out <- list()
  for (tr in table_rows[-1]) {                                  # the first row is the header
    cells <- character(0)
    for (tc in xml2::xml_find_all(tr, "./w:tc")) cells <- c(cells, cell_text(tc))
    if (length(cells) < 11 || !startsWith(cells[1], "rs")) next # empty or auxiliary rows
    location <- regmatches(cells[2], regexec("^(chr[0-9XY]+):([0-9]+)$", cells[2]))[[1]]
    if (length(location) == 0) stop(cells[1], ": could not read the position '", cells[2], "'")
    if (!grepl("^[ACGT]$", cells[4])) stop(cells[1], ": risk allele '", cells[4], "' is not one letter")
    out[[length(out) + 1]] <- c(cells[1], location[2], location[3], cells[3], cells[4], cells[5:11])
  }
  snp_table <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
  names(snp_table) <- c("snp", "chr", "pos", "gene", "risk_allele", "trait", "study",
                        "n_genotyped", "genotype_freq", "hwe_p",
                        "n_tested_probes_paper", "n_sig_probes_paper")
  stopifnot(nrow(snp_table) == 52, !anyDuplicated(snp_table$snp))
  write.csv(snp_table, out_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("%d SNPs written to %s\n", nrow(snp_table), out_file))
}
