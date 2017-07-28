#######################################################################
## This code enables users to get PICS credible sets for nearly any SNP with an index pvalue
## Then users can test if the signals driving the two associations colocalize

### Example workflow:
####R> gwas.credSet <- download.pics(rsid=rs123, pvalue=200)
####R> eQTL.credSet <- read.pics("/full/path/to/pics.txt")
####R> myColoc <- coloc.pics(gwas.credSet, eQTL.credSet)

## Note: more info for each function below
warning(paste("If you are behind a firewall and plan to use download.pics() make sure to set RCurlOptions","Example: options(RCurlOptions = list(proxy=\"gskproxy.gsk.com:800\", proxyuserpwd=\"USER:PW\", useragent = \"R-coloc.pics\"))", sep="\n"), call. = FALSE)
require(data.table)
require(RCurl)
#######################################################################
# coloc.pics : Test for colocalization of two PICS sets
## Example: myColoc <- coloc.pics(gwas.credSet, eQTL.credSet)
## Returns full coloc object. 
##      For "lighter" dataframe returing *only* H3 & H4 posteriors use coloc.pics.lite()
## options:
### data1, data2 = PICS sets from read.pics or download.pics
### pics1, pics2 = column name to pull PICS prob from. Default = "PICS_probability"
### rsid1, rsid2 = column name to pull rsid from.      Default = "Linked_SNP"
### rounded  = [#] Decimal points to round posteriors to
### priorc1  = Prior probability for colocalization with siganl for data1  Default = 1e-4
### priorc2  = Prior probability for colocalization with siganl for data2 Default = 1e-4
### priorc12 = Prior probability for colocalization of both signals.   Default = 1e-5
# Note: Adapted from Toby Johnson's code which was adapted from Giambartolomei et al. 2014 colocalisation method
coloc.pics <- function(data1, 
                       data2,
                       pics1 = "PICS_probability",
                       pics2 = "PICS_probability",
                       rsid1 = "Linked_SNP",
                       rsid2 = "Linked_SNP",
                       rounded = 6,
                       priorc1 = 1e-4, 
                       priorc2 = 1e-4, 
                       priorc12 = 1e-5
                       ) {
  stopifnot(exists("data1") & exists("data2"))
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
# read.pics : load a PICs file that was either download from the website (copy+paste to a text file) or a saved output from download.pics(output=X)
## Example: myCredSet <- read.pics("/full/path/to/PICs_credSet.txt")
read.pics <- function(x){
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
# download.pics : Query the PICs website and get a PICs credible set for a SNP of interest
## Example: myCredSet <- download.pics(rsid="rs6795744", pvalue = "20")
### options:
### rsid = SNP rsID 
### pvalue = -log(pvalue). Example GWAS pvalue 1E^-20 means use pvalue=20
### ancestry = [EUR, ASN, AFR] Default = EUR
### output = Full path and name for the download. Default = temporary file.

download.pics <- function(rsid, pvalue, ancestry = "EUR", output = NA, override = FALSE){
  stopifnot(exists("rsid") & exists("pvalue"))
  if(grepl("\\d+e\\-", "12", ignore.case =TRUE)){
    sub("\\d+e\\-(\\d+)$", "\\1", pvalue, ignore.case =TRUE)
  }
  picsFile <- if(!is.na(output)) output else tempfile() 
  f <- CFILE(picsFile, mode="wb")
  url <- paste('http://pubs.broadinstitute.org/cgi-bin/finemapping/picscgi.pl?command1=', rsid, '&command2=', pvalue, '&command3=', ancestry, sep = "")
  curlPerform(url = url, writedata = f@ref)
  close(f)
  res <- read.pics(picsFile)
  if(!length(res$Linked_SNP)>1 & !override){
    warning("PICs returned only 1 causal SNP. This is most likely b/c the query rsID wasn't in 1KGp1 used for PICs.\n To confirm that there is 1 causal SNP use the PICs website and HaploReg.\n If you are sure 1 causal SNP is correct, use override=TRUE")
    return(NULL)
  }
  if (is.na(output)) unlink(picsFile)
  return(res)
}
#######################################################################

#######################################################################
# coloc.pics.lite : Test for colocalization of two PICS sets
## Example: myColoc <- coloc.pics.lite(gwas.credSet, eQTL.credSet)
## Returns only H3 & H4 posteriors
## options:
### data1, data2 = PICS sets from read.pics or download.pics
### pics1, pics2 = column name to pull PICS prob from. Default = "PICS_probability"
### rsid1, rsid2 = column name to pull rsid from.      Default = "Linked_SNP"
### rounded  = [#] Decimal points to round posteriors to
### priorc1  = Prior probability for colocalization with siganl for data1  Default = 1e-4
### priorc2  = Prior probability for colocalization with siganl for data2 Default = 1e-4
### priorc12 = Prior probability for colocalization of both signals.   Default = 1e-5
# Note: Adapted from Toby Johnson's code which was adapted from Giambartolomei et al. 2014 colocalisation method
coloc.pics.lite <- function(data1, 
                       data2,
                       pics1 = "PICS_probability",
                       pics2 = "PICS_probability",
                       rsid1 = "Linked_SNP",
                       rsid2 = "Linked_SNP",
                       rounded = 6,
                       priorc1 = 1e-4, 
                       priorc2 = 1e-4, 
                       priorc12 = 1e-5
) {
  stopifnot(exists("data1") & exists("data2"))
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
    tmp[[pp1]][n] <- if(!is.na(match(ids[n], data1[[opts$rsid1]]))) data1[which(data1[[opts$rsid1]]==ids[n]),][[opts$pics1]] else 0
    tmp[[pp2]][n] <- if(!is.na(match(ids[n], data2[[opts$rsid1]]))) data2[which(data2[[opts$rsid2]]==ids[n]),][[opts$pics2]] else 0 
  }
  res <- as.data.frame(cbind(
    .norm1(tmp[[pp1]], log = FALSE),
    .norm1(tmp[[pp2]], log = FALSE)
  ))
  colnames(res) <- c(pp1, pp2)
  rownames(res) <- ids
  return(res)
}
