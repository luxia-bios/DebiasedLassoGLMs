# DebiasedLassoGLMs

This repository documents the R code accompanying the article "De-biased Lasso for Generalized Linear Models with A Diverging Number of Covariates" by Lu Xia, Bin Nan and Yi Li (2021+). 

Citation: Xia, L., Nan, B. and Li, Y. (2021+). De-biased Lasso for Generalized Linear Models with A Diverging Number of Covariates. Biometrics, in press.

The R code file `DBL_GLMs_functions.R` provides two functions for inference in generalized linear models (GLMs), one for implementing the proposed de-biased lasso approach by directly inverting the Hessian matrix, and the other for implementing the original de-biased lasso approach (van de Geer et al., 2014). The R code file `DBL_GLMs_example.R` provides some example code to compare the aforementioned two de-biased lasso methods and the maximum likelihood estimation (MLE) with simulated data. For more detailed description, please see the document `README.pdf`.


