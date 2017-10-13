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
    help = "Comma separated: rsID,pval,ethnicity\n\t\tNote: *Header is mandatory & case sensitive.*\n\t\tNote: Header must contain rsID and pval, ethnicity is optional, unlimited other fields are optional."
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
suppressMessages(library(RSQLite, lib.loc="/home/kbs14104/R"))
suppressMessages(library(BiocGenerics, lib.loc="/home/kbs14104/R"))
suppressMessages(library(S4Vectors, lib.loc="/home/kbs14104/R"))
suppressMessages(library(IRanges, lib.loc="/home/kbs14104/R"))
suppressMessages(library(Biobase, lib.loc="/home/kbs14104/R"))
suppressMessages(library(AnnotationDbi, lib.loc="/home/kbs14104/R"))
suppressMessages(library(GenomeInfoDb, lib.loc="/home/kbs14104/R"))
suppressMessages(library(GenomicRanges, lib.loc="/home/kbs14104/R"))
suppressMessages(library(GenomicFeatures, lib.loc="/home/kbs14104/R"))
suppressMessages(library(OrganismDbi, lib.loc="/home/kbs14104/R"))
suppressMessages(library(Homo.sapiens, lib.loc="/home/kbs14104/R"))
suppressMessages(library(dplyr, lib.loc="/home/kbs14104/R"))
suppressMessages(library(dtplyr, lib.loc="/home/kbs14104/R"))
suppressMessages(library(biomaRt, lib.loc="/home/kbs14104/R"))
suppressMessages(source("/home/kbs14104/scripts/coloc_pics/coloc_pics.R"))

ensembl = useMart("ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", dataset = "hsapiens_gene_ensembl")
snp_mart = useMart("ENSEMBL_MART_SNP", host = "grch37.ensembl.org", dataset="hsapiens_snp")
snp_attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end" )

eqtl.files <- read.table(file="/GWD/appbase/projects/statgen3/PICS_sets/ToolVariant_cred-sets/all_eQTLs_filePaths_table.csv", 
                         header = TRUE, 
                         sep = ",", 
                         stringsAsFactors = FALSE)
eqtl.files <- as.data.table(eqtl.files)
setkey(eqtl.files, geneID)
eqtl.tissues <- unique(eqtl.files[["tissue"]])

# Known GWAS loci to look up
known.Loci <- read.table(file=opt$input, sep = ",", header = TRUE, stringsAsFactors=FALSE, fileEncoding="latin1")

# DF to keep final results in
known.Loci.coloc <- as.data.frame(matrix(ncol = length(colnames(known.Loci))+6))
# testing new Subset here
known.Loci.coloc <- as.data.frame(matrix(ncol = eval(length(colnames(known.Loci))+6), nrow =  eval(length(eqtl.tissues) * dim(known.Loci)[1] * 50)))
colnames(known.Loci.coloc) <- c(colnames(known.Loci), "tissue", "eqtl_index_snp", "hgnc_symbol", "ensembl_ID", "H3", "H4")
# row counter
r <- 1 

for(i in 1:dim(known.Loci)[1]){
  cat("\tColocalizing ", known.Loci[i,]$rsID, " at ", known.Loci[i,]$pval, sep = "", file=stderr())
  # Determine ancestry for PICs
  ethnicity <- "EUR"
  if( ("ethnicity" %in% colnames(known.Loci)) ){
    if (grepl("EAS", known.Loci[i,]$ethnicity)){
      ethnicity <- "ASN"
    }
    else if (grepl("AA", known.Loci[i,]$ethnicity)){
      ethnicity <- "AFR"
    }
  }
  # Get & clean PICs from website for GWAS
  gwas.pics <- download.pics(rsid = known.Loci[i,]$rsID, pvalue = known.Loci[i,]$pval, ancestry = ethnicity)
  if(is.null(gwas.pics)){next}
  gwas.pics <- gwas.pics[grep("^rs", gwas.pics$Linked_SNP, perl = TRUE),]
  gwas.bed <- getBM(values=gwas.pics$Linked_SNP, 
                    attributes=snp_attributes, 
                    filters="snp_filter", 
                    mart=snp_mart)
  colnames(gwas.bed)  <- c("rsID", "chr", "start", "end")
  gwas.bed$chr <- paste("chr", gwas.bed$chr, sep="")
  # Determine the region for +/- 500 kb
  credSet.limits <- data.frame(chrom=gwas.bed[1,]$chr, 
                               start=(min(gwas.bed$start)-500000), 
                               end=(max(gwas.bed$end)+500000)) %>% makeGRangesFromDataFrame
  # Get all genes in +/- 500 kb
  credSet.genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), credSet.limits) %>% as.data.frame
  if(!dim(credSet.genes)[1]>0){next}
  # Append gene names
  genemap <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"), 
                   filters = "entrezgene",
                   values = credSet.genes$gene_id, 
                   mart = ensembl 
  )
  idx <- match(credSet.genes$gene_id, genemap$entrezgene)
  credSet.genes$hgnc_symbol <- genemap$hgnc_symbol[idx]
  credSet.genes$ensembl_gene_id <- genemap$ensembl_gene_id[idx]
  credSet.genes <- credSet.genes[!is.na(credSet.genes$hgnc_symbol),]
  if(!length(credSet.genes$hgnc_symbol)>0){next}
  
  locus.eqtl.info <- na.omit(eqtl.files[.(credSet.genes$ensembl_gene_id)], cols="filePath")
  if(!dim(locus.eqtl.info)[1]>=1){next}

  z <- 1
  for (e in 1:dim(locus.eqtl.info)[1]){
    if(z==10){
      cat(" .", file=stderr(), sep = "")
      z <- 1
    }
    # Skip if the eQTL cred set file is missing
    if(!file.exists(locus.eqtl.info[e]$filePath)){next}
    # Load eQTL cred set
    eqtl.pics <- read.pics(locus.eqtl.info[e]$filePath)
    # Skip if empty
    if(dim(eqtl.pics)[1]<=1){next}
    # Calc and append/return GWAS:eQTL coloc
    known.Loci.coloc[r,] <- c(known.Loci[i,], 
        locus.eqtl.info[e]$tissue, 
        locus.eqtl.info[e]$index.rsID, 
        locus.eqtl.info[e]$hgnc_symbol, 
        locus.eqtl.info[e]$geneID, 
        coloc.pics.lite(gwas.pics, eqtl.pics))
    r <- r+1
    z <- z+1
  }
  cat("\n", sep="", file=stderr())
}


known.Loci.coloc <- known.Loci.coloc[!is.na(known.Loci.coloc$rsID),]

write.csv(known.Loci.coloc, file=opt$output, row.names=FALSE, quote=FALSE)
cat("FINISHED coloc analysis\n", sep="", file=stderr())