#!/home/kbs14104/bin/Rscript
readRenviron("~/.Renviron")
.libPaths(c("/home/kbs14104/R/", 
  "/home/kbs14104/R/x86_64-pc-linux-gnu-library/3.3", 
  "/home/kbs14104/etc/R-3.3.3/library", 
  .libPaths())
)
suppressMessages(library(optparse, lib.loc="/home/kbs14104/R"))
options(warn=0)

option_list = list(
  make_option(
    c("-i", "--input"),
    default = NULL,
    metavar = "FILE.csv",
    help = "Comma separated: rsID,pval,ancestry\n\t\tNote: *Header is mandatory & case sensitive.*\n\t\tNote: Header must contain rsID and pval, ancestry is optional, unlimited other fields are optional."
  ),
  make_option(
    c("-o", "--output"),
    default = NULL,
    metavar = "FILE.csv",
    help = "</path/for/output/example.csv>"
  )
)
opt_parser = optparse::OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Must use --input=example.csv\n", call. = FALSE)
}
if(!file.exists(opt$input)){
  stop("Input does not exist.\n", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Must use --output=example.csv\n", call. = FALSE)
}

suppressMessages(library(data.table, lib.loc="/home/kbs14104/R"))
suppressMessages(library(plyr,       lib.loc="/home/kbs14104/R"))
suppressMessages(library(tidyverse,  lib.loc="/home/kbs14104/R"))
suppressMessages(library(stringr,    lib.loc="/home/kbs14104/R"))
suppressMessages(library(RSQLite,    lib.loc="/home/kbs14104/R"))
# suppressMessages(source("/GWD/appbase/projects/RD-Genetic-INSiGHT/PICCOLO/scripts/coloc_pics/coloc_pics.R"))
suppressMessages(source("/home/kbs14104/scripts/coloc_pics/coloc_pics.R"))

piccolo_db <- DBI::dbConnect(
  RSQLite::SQLite(), 
  "/GWD/appbase/projects/RD-Genetic-INSiGHT/PICCOLO/PICCOLO_test.sqlite")

snp150_tbl        <- tbl(piccolo_db, "snp150_hg19")
gtexGene_tbl      <- tbl(piccolo_db, "gtexGene")
piccolo_eqtls_tbl <- tbl(piccolo_db, "piccolo_eqtl_files")

# Known GWAS loci to look up
input <- read_csv(file=opt$input, col_types = cols(), col_names = TRUE, progress = FALSE)
# Check input columns
if(is.null(input[["rsID"]]))     { stop("Input is missing column: rsID\n", call. = FALSE) }
if(is.null(input[["pval"]]))     { stop("Input is missing column: pval\n", call. = FALSE) }
if(is.null(input[["ancestry"]])) { input$ancestry <- "EUR" }

# clean input names to avoid join by clash with other 
input <- input %>% 
  plyr::rename(replace = c("chrom" = "chrom_input"), warn_missing = FALSE) %>% 
  plyr::rename(replace = c("chr"   = "chr_input"),   warn_missing = FALSE) %>% 
  plyr::rename(replace = c("pos"   = "pos_input"),   warn_missing = FALSE)

ret <- input %>% 
  mutate(rsID_short = str_replace(rsID, "rs(\\d+)", "\\1") %>% as.integer) 

## First, get hg19 coords for input snps
gwas <- inner_join(snp150_tbl, 
             ret, 
             by = c("id" = "rsID_short"), 
             copy = TRUE) %>% collect() %>%   
  select(-rsid) %>% 
  ## Next, filter for genes with 500 kb
  group_by_all() %>% 
  do(tbl(piccolo_db, sql(paste0("SELECT * FROM gtexGene WHERE chr == ", .$chr, " AND (",
        "(start >= ",     .$start - 5e5, " AND start <= ", .$end + 5e5, ")",
        " OR (end >= ",   .$start - 5e5, " AND end <= ",   .$end + 5e5, ")", 
        " OR (start >= ", .$start - 5e5, " AND end <= ",   .$end + 5e5, ")",
        ")"))) %>% 
      collect() %>% 
      select(hgnc_id, ensembl_id)) %>% 
  ## For the nearby genes, find PICs cred set files
  inner_join(piccolo_eqtls_tbl, ., by = c("ensembl_idx" = "ensembl_id"), copy = TRUE, suffix = c("", "_input")) %>% 
  collect()

## Download PICs for input GWAS 
gwas <- inner_join(gwas, 
                   input %>% 
                     group_by_at(.vars = names(input)) %>% 
                     do(pics_gwas = pics.download(rsid     = .$rsID, 
                                                  pvalue   = .$pval, 
                                                  ancestry = .$ancestry)),
                   by = names(input)) %>% 
        select(-chr, -chrom, -start, -end, -rsid_idx_short)

## Run colocs
gwas <- gwas %>% 
  group_by_at(.vars = c(names(input), "rsid_idx", "tissue_idx", "hgnc_idx", "ensembl_idx", "pval_idx", "pmid_idx", "filePath_idx")) %>% 
  do(pics.coloc.lite(.$pics_gwas %>% .[[1]], pics.read(.$filePath_idx))) %>% 
  ungroup()

# Append results with no chrom:pos as NAs
gwas <- right_join(gwas, input, by = names(input)) %>% 
  select(-matches("filePath")) %>% 
  plyr::rename(replace = c("chrom_input" = "chrom"), warn_missing = FALSE) %>% 
  plyr::rename(replace = c("chr_input"   = "chr"),   warn_missing = FALSE) %>% 
  plyr::rename(replace = c("pos_input"   = "pos"),   warn_missing = FALSE)

write_csv(gwas, path = opt$output, col_names = TRUE)

dbDisconnect(piccolo_db)
cat("FINISHED coloc analysis\n", sep="", file=stderr())
