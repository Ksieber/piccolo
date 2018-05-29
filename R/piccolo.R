suppressMessages(library(data.table))
suppressMessages(library(RCurl))
#' pics.coloc : Test for colocalization of two PICS sets
#' Adapted from [Toby Johnson's code](https://github.com/tobyjohnson/gtx/blob/master/R/abf.R "Toby's coloc in R")
#' Which was adapted from [Giambartolomei et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24830394 "Giambartolomei et al. 2014") colocalization method
#' @examples 
#' myColoc <- pics.coloc(gwas.credSet, eQTL.credSet)
#' 
#' @return Full colocalization info. For concise dataframe returing *only* H3 & H4 posteriors use pics.coloc.lite()
#' @param data1  PICS sets, from pics.read or pics.download
#' @param data2  PICS sets, from pics.read or pics.download
#' @param rounded  [#] Decimal points to round posteriors to
#' @param priorc1  Prior probability for colocalization with siganl for data1   Default = 1e-4
#' @param priorc2  Prior probability for colocalization with siganl for data2   Default = 1e-4
#' @param priorc12 Prior probability for colocalization of both signals.        Default = 1e-5
#' @param pics1  Column name to pull PICS prob from data1 if data not from pics.* Default = "PICS_probability"
#' @param pics2  Column name to pull PICS prob from data2 if data not from pics.* Default = "PICS_probability"
#' @param rsid1  Column name to pull rsid from if data not from pics.*            Default = "Linked_SNP"
#' @param rsid2  Column name to pull rsid from if data not from pics.*            Default = "Linked_SNP"
#' @export 
pics.coloc <- function(data1, 
                       data2,
                       pics1    = "PICS_probability", # column header for poster probabilities in data1
                       pics2    = "PICS_probability", # column header for poster probabilities in data2
                       rsid1    = "Linked_SNP",       # column header for snps in LD in data1
                       rsid2    = "Linked_SNP",       # column header for snps in LD in data2
                       rounded = 6,
                       priorc1 = 1e-4, 
                       priorc2 = 1e-4, 
                       priorc12 = 1e-5
                       ) {
  stopifnot(exists("data1") & exists("data2"))
  if(is.logical(data1)){
    if(is.na(data1)){
      warning("pics.coloc- data1 is NA, skipping coloc.\n")
      return(list(results = NA, nvariants = NA))
    }
  }
  if(is.logical(data2)){
    if(is.na(data2)){
      warning("pics.coloc - data2 is NA, skipping coloc.\n")
      return(list(results = NA, nvariants = NA))
    }
  }
  pics <- .harmonize.pics(data1, 
                          data2, 
                          opts <- data.frame(rsid1 = rsid1, rsid2 = rsid2, pics1 = pics1, pics2 = pics2, stringsAsFactors = FALSE))
  
  nv <- dim(pics)[1]
  res <- data.frame(hypothesis = paste0("H", 3:4),
                    label = c("Two variants separately associated with phenotypes 1 and 2",
                              "One variant associated with phenotypes 1 and 2"),
                    prior = .norm1(c(priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf    = c((sum(pics[[1]])*sum(pics[[2]]) - sum(pics[[1]]*pics[[2]]))/(nv*(nv - 1)), 
                            sum(pics[[1]]*pics[[2]])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- .norm1(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(list(results = res, nvariants = nv))
}

#######################################################################
#' pics.read : read a PICs file 
#' 
#' Read into memory PICs from pics.download(output="X") or manually saved text from the website (copy+paste to a text file, not suggested)
#' @examples 
#' myCredSet <- pics.read("/full/path/to/PICs_credSet.txt")
#' 
#' @param filepath String with the filepath to the data to read. 
#' @return PICs dataframe
#' @export
pics.read <- function(x){
  stopifnot(exists("x"))
  stopifnot(file.exists(x))
  # Pull the first 3 lines to check if this file was downloaded with wget (has fluff header) or was copy+pasted from gui (no fluff header)
  fh <- file(x, "r")
  head <- readLines(fh, n=3)
  close(fh)
  if(grepl("^\\<pre\\>", head[1], perl = TRUE) & grepl("^Index_SNP", head[1], perl = TRUE)){
    pics <- fread(file=x, header=TRUE, sep="\t", skip = 2)
  }
  else {
    pics <- fread(file=x, header=TRUE, sep="\t")
  }
  return(pics)
}
#######################################################################


#######################################################################
#' pics.download : Query the PICs website and get a PICs credible set for a SNP of interest
#' @examples 
#' myCredSet <- pics.download(rsid="rs6795744", pvalue = "20")
#' 
#' @param rsid valid rsID. e.g. rs6795744  
#' @param pvalue p-value for the rsid. Either scientific or -log10 base. e.g. "1E-20", "1e-20", or "20"
#' @param ancestry = [EUR, ASN, AFR] Default = EUR
#' @param output Specify path and name for the download to be saved to. Default = temporary file.
#' @return PICs credible set
#' @export
pics.download <- function(rsid, pvalue, ancestry = "EUR", output = NA, override = FALSE){
  stopifnot(exists("rsid") & exists("pvalue"))
  if(grepl("\\d\\.?[eE]\\-\\d+", pvalue, perl=TRUE)){
    pvalue <- -log10(pvalue)
  }
  picsFile <- if(!is.na(output)) output else tempfile() 
  fileHandle <- CFILE(picsFile, mode = "wb")
  url <- paste('http://pubs.broadinstitute.org/cgi-bin/finemapping/picscgi.pl?command1=', rsid, '&command2=', pvalue, '&command3=', ancestry, sep = "")
  curlPerform(url = url, writedata = fileHandle@ref)
  close(fileHandle)
  res <- pics.read(picsFile)
  if(!length(res$Linked_SNP) > 1 & !override){
    if(options("warn") > 0){
      warning("PICs returned only 1 causal SNP. This is most likely b/c the query rsID wasn't in 1KGp1 used for PICs.\n To confirm that there is 1 causal SNP use the PICs website and HaploReg.\n If you are sure 1 causal SNP is correct, use override=TRUE\n");
    }
    return(NA)
  }
  if (is.na(output)) unlink(picsFile)
  return(res)
}
#######################################################################

#######################################################################
#' pics.coloc.lite : Test for colocalization of two PICS sets
#' 
#' Same arguements as \code{\link{pics.coloc}}, return object is the only difference.
#' @examples 
#' myColoc  <- pics.coloc(gwas.credSet, eQTL.credSet)
#' 
#' @examples 
#' myColocs <- data %>% group_by() %>% do(pics.coloc.lite(.$gwas_pics %>% .[[1]], .$eqtl_pics %>% .[[1]]))
#' 
#' @return H3 & H4 posterior colocalization info. For concise dataframe returing *only* H3 & H4 posteriors use pics.coloc.lite()
#' @param data1  PICS sets, from pics.read or pics.download
#' @param data2  PICS sets, from pics.read or pics.download
#' @return H3 & H4
#' @export
pics.coloc.lite <- function(data1, 
                       data2,
                       pics1    = "PICS_probability", # column header for poster probabilities in data1
                       pics2    = "PICS_probability", # column header for poster probabilities in data2
                       rsid1    = "Linked_SNP",       # column header for snps in LD in data1
                       rsid2    = "Linked_SNP",       # column header for snps in LD in data2
                       rounded  = 6,
                       priorc1  = 1e-4, 
                       priorc2  = 1e-4, 
                       priorc12 = 1e-5
) {
  stopifnot(exists("data1") & exists("data2"))
  if(is.logical(data1)){
    if(is.na(data1)){
      if(options("warn")>0){warning("pics.coloc - data1 is NA, skipping coloc.\n");}
      return(data.frame(H3 = NA, H4 = NA))
    }
  }
  if(is.logical(data2)){
    if(is.na(data2)){
      if(options("warn")>0){warning("pics.coloc - data2 is NA, skipping coloc.\n");}
      return(data.frame(H3 = NA, H4 = NA))
    }
  }
  # Progress bar
  # if ((!is.null(.pb)) && inherits(.pb, "Progress") && (.pb$i < .pb$n)) .pb$tick()$print()

  pics <- .harmonize.pics(data1, 
                          data2, 
                          opts <- data.frame(rsid1 = rsid1, rsid2 = rsid2, pics1 = pics1, pics2 = pics2, stringsAsFactors = FALSE))
  
  nv <- dim(pics)[1]
  res <- data.frame(prior = .norm1(c(priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf    = c((sum(pics[[1]])*sum(pics[[2]]) - sum(pics[[1]]*pics[[2]]))/(nv*(nv - 1)), 
                              sum(pics[[1]]*pics[[2]])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- .norm1(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(data.frame(H3=res$posterior[1], H4=res$posterior[2]))
}

#######################################################################
## Internal functions
#######################################################################
.norm1 <- function(x, log = FALSE) {
  if (all(is.na(x))) return(x)
  if (log == TRUE) {
    x <- x - max(x, na.rm = TRUE)
    x <- exp(x)    
  } else {
    ## This does not work if x contains NaNs or +Infs
    stopifnot(all(x >= 0, na.rm = TRUE))
    x <- x / max(x, na.rm = TRUE)
  }
  return(x / sum(x, na.rm = TRUE))
}

.harmonize.pics <- function(data1, 
                            data2, 
                            opts = data.frame(pics1 = "PICS_probability",
                                              pics2 = "PICS_probability",
                                              rsid1 = "Linked_SNP",
                                              rsid2 = "Linked_SNP",
                                              stringsAsFactors = FALSE)
                            ){
  ids <- unique(c(data1[[opts$rsid1]], data2[[opts$rsid2]]))
  tmp <- as.data.frame(matrix(data = NA, nrow = length(ids), ncol = 2))
  pp1 <- if (opts$pics1==opts$pics2) paste(opts$pics1, ".1", sep = "") else opts$pics1
  pp2 <- if (opts$pics1==opts$pics2) paste(opts$pics2, ".2", sep = "") else opts$pics2
  colnames(tmp) <- c(pp1, pp2)
  for(n in 1:length(ids)){
    tmp[[pp1]][n] <- if(!is.na(match(ids[n], data1[[opts$rsid1]]))) data1[which(data1[[opts$rsid1]]==ids[n]),][[opts$pics1]][1] else 0
    tmp[[pp2]][n] <- if(!is.na(match(ids[n], data2[[opts$rsid1]]))) data2[which(data2[[opts$rsid2]]==ids[n]),][[opts$pics2]][1] else 0 
  }
  res <- as.data.frame(cbind(
    .norm1(tmp[[pp1]], log = FALSE),
    .norm1(tmp[[pp2]], log = FALSE)
  ))
  colnames(res) <- c(pp1, pp2)
  rownames(res) <- ids
  return(res)
}
