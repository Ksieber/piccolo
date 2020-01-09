# piccolo
This code enables users to get PICS credible sets for nearly any index SNP with a pvalue. 
Users can then test if two genetic signals colocalize in a Bayesian framework. 
***

## Example workflow:
```
library("piccolo")
gwas_credSet <- pics.download(rsid="rs123", pvalue=-log10(1.23E-45))
eQTL_credSet <- pics.read("/full/path/to/pics.txt")
myColoc      <- pics.coloc(gwas_credSet, eQTL_credSet)
```
See the pinned issue for a pptx showing further illustrating an example.  Also available [here](https://github.com/Ksieber/piccolo/files/4041189/piccolo_info.pptx "Example slide deck")
***


## pics.coloc: 
Test for colocalization of two PICS credible sets
#### Example: `myColoc <- pics.coloc(gwas_credSet, eQTL_credSet)`   
#### Options  
|Options      | Value 						   | Default
| ----------- | ---------------------------------------------------------- | ------------------ |
|data1, data2 | PICS sets from pics.read or pics.download 				   | NA
|pics1, pics2 | column name to pull PICS prob from. 					   | "PICS_probability"
|rsid1, rsid2 | column name to pull rsid from. 							   | "Linked_SNP"
|rounded  	  | [#] Decimal points to round posteriors to. 				   |
|priorc1  	  | Prior probability for a siganl in data1 | 1e-4
|priorc2  	  | Prior probability for a siganl in data2 | 1e-4 
|priorc12 	  | Prior probability for colocalization of both signals.      | 1e-5   
#### Note:   
 * Returns full coloc object.  
 * Users only need to specify the column names through pics[1|2] & rsid[1|2] if they are using credible sets that are *not* from pics.download().  
 * For "lighter" dataframe returing *only* H3 & H4 posteriors, use pics.coloc.lite() 
 * Adapted from [Toby Johnson's code](https://github.com/tobyjohnson/gtx/blob/master/R/abf.R "Toby's coloc in R"), which was adapted from [Giambartolomei et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24830394 "Giambartolomei et al. 2014") colocalisation method  
***


## pics.coloc.lite: 
Test for colocalization of two PICS credible sets
#### Example: `myColoc <- pics.coloc.lite(gwas_credSet, eQTL_credSet)`
#### Options: Same as pics.coloc()
#### Note:  
 * Returns a "lighter" dataframe with *only* H3 & H4 posteriors  
***


## pics.download:  
Query the [PICs website](http://pubs.broadinstitute.org/pubs/finemapping/ "PICs") and get a PICs credible set for a SNP of interest
#### Example: `myCredSet <- pics.download(rsid="rs6795744", pvalue = "1.23e-20")`
#### Options:
|Options	| Value | Default
| --------- | ------------------------------------------------ | ----------------------------------- |
|rsid 		| SNP rsID 						  				   | Example: "rs1234"
|pvalue 	| pvalue or -log10(pvalue). 				   | Examples: -log10(1.23E-45) or 45
|ancestry 	| [EUR, ASN, AFR] 				  				   | "EUR"
|output 	| Full path and name for the download to be saved. | held in memory
#### Note:
 * Returns a dataframe with default colnames as pic[1|2]="PICS_probability" & rsid[1|2]="Linked_SNP" that pics.coloc anticipates as default. 
 * For pics.download() to work behind a proxy use this code in an R session in your ~/.Renviron. 
 ```
 options(RCurlOptions = list(proxy="proxy.URL:PORT", proxyuserpwd="ID:PW"))
 ```
***

## pics.read: 
Load a PICs file that was either download from the website (copy+paste to a text file) or a saved output from pics.download(output=X)
#### Example: `myCredSet <- pics.read("/full/path/to/PICs_credSet.txt")`
***

## Known Bugs:  
 * If the rsID is not in the imputation used by PICs (1KGp1) then it returns the input rsID as the 1 causal SNP. This is almost certainly an erroneous credible set. This will also occur with rsIDs that do not exist. By default, pics.download() return NA.   
