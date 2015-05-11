# compgen
R files for Computational Genomics final project

All code was written in R.  The following R files are available on GitHub:

- project_code_MatrixFact: This file contains all code used to perform the matrix factorization method and evaluate on simulated and actual data. Code used to generate the simulated data is also provide, but commented out, as the data set is provided.  Code used to evaluate the method on actual data is also commented out.  Running this code with the data_MatrixFact.RData file will test the method on the simulated data.

- data_MatrixFact: This data file contains simulated data (samplesSim) and the random matrix used to generate the data (matrix100).  These data sets can be used to run the project_code_MatrixFact.R file.

- project_code_PCA: This file contains all code used to perform the PCA method and to evaluated on simulated data.  Code used to generate the simulated data is also provided, but commented out, as the data set is provided.  Code used to evaluated this method with three cell types (two dimensions) is also commented out.  Running this code with the data_PCA.RData will test the method on the simulated data with 6 cell types.

- data_PCA: This data file contains simulated data (samples500) and the random matrix used to generate the data (matrix500).  These data sets can be used to run the projet_code_PCA.R file.
