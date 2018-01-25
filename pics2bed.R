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
    metavar = "PICs.txt",
    help = "Saved output *from* pull_PICs.pl"
  ),
  make_option(
    c("-o", "--output"),
    default = NULL,
    metavar = "PICs.bed",
    help = "</path/for/output/example.bed>"
  ),
  make_option(
    c("-r", "--rsid"),
    default = NULL,
    metavar = "rs12345",
    help = "Use rsid and pval instead of saved PICs.txt file"
  ),
  make_option(
    c("-p", "--pval"),
    default = NULL,
    type = "numeric",
    metavar = "1.23E-45",
    help = "Use rsid and pval instead of saved PICs.txt file"
  ),
  make_option(
    c("-e", "--ethnicity"),
    default = "EUR",
    metavar = "ethnicity",
    help = "EUR|ASN|AFR"
  ),
    make_option(
    c("-c", "--cutoff"),
    default = "1",
    type = "numeric",
    metavar = "(0:1]",
    help = "PICs cummalative probability"
  ),
  make_option(
    "--host",
    default = "grch37.ensembl.org",
    metavar = "grch37.ensembl.org",
    help = "hg19=grch37.ensembl.org or hg38=ensembl.org"
  )
)
opt_parser = optparse::OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) & is.null(opt$rsid)) {
  print_help(opt_parser)
  stop("Must use --input=PICs.txt\n", call. = FALSE)
}

if(!is.null(opt$input)){
  if(!file.exists(opt$input)){
    stop("Input does not exist.\n", call. = FALSE)
  }
} 

if(!is.null(opt$rsid) & is.null(opt$pval)){
  stop("--pval must be defined with --rsid.\n", call. = FALSE)
}

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Must use --output=PICs.bed\n", call. = FALSE)
}

suppressMessages(library(biomaRt))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(source("/home/kbs14104/scripts/coloc_pics/coloc_pics.R"))

snp_mart <- useMart("ENSEMBL_MART_SNP", host = opt$host, dataset="hsapiens_snp")
snp_attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end" )

if(!is.null(opt$input)){
  cat("\tReading in PICs credible set\n", sep="", file=stderr())
  picsIn <- read.pics(opt$input)
} else {
  cat("\tCalculating PICs credible set\n", sep="", file=stderr())
  picsIn <- download.pics(rsid = opt$rsid, pvalue = -log10(opt$pval), ancestry = opt$ethnicity)
}

picsIn <- filter(picsIn, cumsum(picsIn$PICS_probability)<=opt$cutoff)

cat("\tConverting input to bed\n", sep="", file=stderr())
bedOut <- getBM(values=picsIn$Linked_SNP, 
      attributes=snp_attributes, 
      filters="snp_filter", 
      mart=snp_mart) %>%
  mutate(chr_name = sub("^(\\d+)$", "chr\\1", .$chr_name, perl=TRUE)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE, 
                           seqnames.field = "chr_name", 
                           start.field = "chrom_start", 
                           end.field = "chrom_end") %>% 
  sort.GenomicRanges %>%
  as.data.frame %>%
  mutate(start=.$start-1) %>%
  subset(select = c(-width, -strand))

bedOut$PICs <- picsIn[match(bedOut$refsnp_id, picsIn$Linked_SNP),]$PICS_probability

write.table(bedOut, file=opt$output, sep="\t", row.names=F, col.names=F, quote=F)
cat("\tFINISHED pics2bed.R\n", sep="", file=stderr())