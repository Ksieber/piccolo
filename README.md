# coloc_pics
This code enables users to get PICS credible sets for nearly any index SNP with a pvalue. 
Then users can test if two genetic signals colocalize in a Bayesian framework. 
***

## Example workflow:
```
source("coloc_pics.R")
gwas.credSet <- download.pics(rsid="rs123", pvalue="1.23E-45")
eQTL.credSet <- read.pics("/full/path/to/pics.txt")
myColoc <- coloc.pics(gwas.credSet, eQTL.credSet)
```
***


## coloc.pics: 
Test for colocalization of two PICS credible sets
#### Example: `myColoc <- coloc.pics(gwas.credSet, eQTL.credSet)`   
#### Options  
|Options      | Value 						   | Default
| ----------- | ---------------------------------------------------------- | ------------------ |
|data1, data2 | PICS sets from read.pics or download.pics 				   | NA
|pics1, pics2 | column name to pull PICS prob from. 					   | "PICS_probability"
|rsid1, rsid2 | column name to pull rsid from. 							   | "Linked_SNP"
|rounded  	  | [#] Decimal points to round posteriors to. 				   |
|priorc1  	  | Prior probability for colocalization with siganl for data1 | 1e-4
|priorc2  	  | Prior probability for colocalization with siganl for data2 | 1e-4 
|priorc12 	  | Prior probability for colocalization of both signals.      | 1e-5   
#### Note:   
 * Returns full coloc object.  
 * Users only need to specify pics[1|2] & rsid[1|2] if they are using credible sets that are *not* from download.pics().  
 * For "lighter" dataframe returing *only* H3 & H4 posteriors use coloc.pics.lite() 
 * Adapted from [Toby Johnson's code](https://github.com/tobyjohnson/gtx/blob/master/R/abf.R "Toby's coloc in R"), which was adapted from [Giambartolomei et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24830394 "Giambartolomei et al. 2014") colocalisation method  
***


## coloc.pics.lite: 
Test for colocalization of two PICS credible sets
#### Example: `myColoc <- coloc.pics.lite(gwas.credSet, eQTL.credSet)`
#### Options: Same as coloc.pics()
#### Note:  
 * Returns a "lighter" dataframe with *only* H3 & H4 posteriors  
***


## download.pics:  
Query the [PICs website](http://pubs.broadinstitute.org/pubs/finemapping/ "PICs") and get a PICs credible set for a SNP of interest
#### Example: `myCredSet <- download.pics(rsid="rs6795744", pvalue = "1.23e-20")`
#### Options:
|Options	| Value | Default
| --------- | ------------------------------------------------ | ----------------------------------- |
|rsid 		| SNP rsID 						  				   | Example: rs1234
|pvalue 	| -log(pvalue) or full GWAS pval. 				   | Examples: "1E-20", "1e-20", "20"
|ancestry 	| [EUR, ASN, AFR] 				  				   | EUR
|output 	| Full path and name for the download to be saved. | held in memory
|override | If TRUE, use cred sets with 1 causal SNP. *CAUTION* using this. Read below "Known Bugs". | FALSE
#### Note:
 * Returns a dataframe with default colnames as pic[1|2]="PICS_probability" & rsid[1|2]="Linked_SNP" that coloc.pics anticipates as default. 
 * For download.pics() to work behind a proxy use this code in an R session in your ~/.Renviron. 
 ```
 options(RCurlOptions = list(proxy="proxy.URL:PORT", proxyuserpwd="ID:PW", useragent = "R-coloc.pics"))
 ```
***

## read.pics : 
Load a PICs file that was either download from the website (copy+paste to a text file) or a saved output from download.pics(output=X)
#### Example: `myCredSet <- read.pics("/full/path/to/PICs_credSet.txt")`
***

## Known Bugs:  
 * If the rsID is not in the imputation used by PICs (1KGp1) then it returns the input rsID as the 1 causal SNP. This is almost certainly an erroneous credible set. This will also occur with rsIDs that do not exist. By default, download.pics() will skip *any* cred set with 1 causal SNP to avoid problems with missing/fake data. If you are certain that 1 causal SNP is correct for *a* locus, set override=TRUE in download.pics().  
 
 ## colocPics_GTEx.R
 This script will take input genetic signals and colocalize them with all GTEx eQTLs. This can only be run on the GSK filesystem currently. Use --help for options. 
