# coloc_pics
This code enables users to get PICS credible sets for nearly any index SNP with a pvalue. 
Then users can test if two genetic signals colocalize in a Bayesian framework. 

***
## Example workflow:
```
gwas.credSet <- download.pics(rsid=rs123, pvalue=200)
eQTL.credSet <- read.pics("/full/path/to/pics.txt")
myColoc <- coloc.pics(gwas.credSet, eQTL.credSet)
```
***


***
## coloc.pics : Test for colocalization of two PICS credible sets
Example: `myColoc <- coloc.pics(gwas.credSet, eQTL.credSet)`
*Returns full coloc object. 
*For "lighter" dataframe returing *only* H3 & H4 posteriors use coloc.pics.lite()
####Options
|Options	  | Value 													   | Default
| ----------- | ---------------------------------------------------------- | ------------------ |
|data1, data2 | PICS sets from read.pics or download.pics 				   | NA
|pics1, pics2 | column name to pull PICS prob from. 					   | "PICS_probability"
|rsid1, rsid2 | column name to pull rsid from. 							   | "Linked_SNP"
|rounded  	  | [#] Decimal points to round posteriors to. 				   |
|priorc1  	  | Prior probability for colocalization with siganl for data1 | 1e-4
|priorc2  	  | Prior probability for colocalization with siganl for data2 | 1e-4
|priorc12 	  | Prior probability for colocalization of both signals.      | 1e-5
####Note: Adapted from [Toby Johnson's code](https://github.com/tobyjohnson/gtx/blob/master/R/abf.R "Toby's coloc in R"), which was adapted from [Giambartolomei et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24830394 "Giambartolomei et al. 2014") colocalisation method
***


***
## coloc.pics.lite : Test for colocalization of two PICS credible sets
Example: `myColoc <- coloc.pics.lite(gwas.credSet, eQTL.credSet)`
	..*For "lighter" dataframe returing *only* H3 & H4 posteriors
####Options: Same as coloc.pics()
***


***
## download.pics : Query the PICs website and get a PICs credible set for a SNP of interest
Example: `myCredSet <- download.pics(rsid="rs6795744", pvalue = "20")`
####Options:
|Options	| Value | Default
| --------- | ------------------------------------------------ | ----------------------------------- |
|rsid 		| SNP rsID 						  				   | Example: rs1234
|pvalue 	| -log(pvalue) or full GWAS pval. 				   | Examples: "1E-20", "1e-20", "20"
|ancestry 	| [EUR, ASN, AFR] 				  				   | Default = EUR
|output 	| Full path and name for the download to be saved. | Default = temporary file in memory.
***


***
## read.pics : load a PICs file that was either download from the website (copy+paste to a text file) or a saved output from download.pics(output=X)
Example: `myCredSet <- read.pics("/full/path/to/PICs_credSet.txt")`
***