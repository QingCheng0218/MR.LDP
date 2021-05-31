MR.LDP
=======

**MR.LDP** is a package for two-sample Mendelian randomization for GWAS summary statistics accounting linkage disequilibrium and horizontal pleiotropy(MR-LDP).

Installation
============
Install the development version of **MR.LDP** by use of the 'devtools' package. Note that **MR.LDP** depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```
library(devtools)
install_github("QingCheng0218/MR.LDP")
```

If you have errors installing this package on Linux, try the following commands in R:
```
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252") 
Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(devtools)
install_github("QingCheng0218/MR.LDP@main")
```

Usage
=========
The ['MR.LDP' vignette](https://github.com/QingCheng0218/MR.LDP/blob/master/vignettes/MRLDP.pdf) will provide a good start point for two-sample Mendelian randomization analysis using **MR.LDP** package. 

References
==========
Qing Cheng, Yi Yang, Xingjie Shi, Kar-Fu Yeung, Can Yang, Heng Peng, Jin Liu<sup>+</sup>. [MR-LDP: a two-sample Mendelian randomization for GWAS summary statistics accounting linkage disequilibrium and horizontal pleiotropy. NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020.](https://academic.oup.com/nargab/article/2/2/lqaa028/5828855?login=true)

Development
===========

This package is developed and maintained by Qing Cheng (qing.cheng@duke-nus.edu.sg).
 
