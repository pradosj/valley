valley
======
valley is an R package containing various functions for signal processing.
The package was originaly tested for peak detection in mass spectrometry (1D, 2D, 3D signals), but the functions 
are not limited to this domain, and can for example be useful for spot detection in images (see package example).

The core part of the algorithms is implemented in c++11 with Rcpp framework. The package those depends on Rcpp package, 
and also makes use of roxygene2 to generate its documentation from source code.


Installation
---------------
To install the package you can either load the Rstudio project and run Build, or run the commands:

    Rcpp::compileAttributes()
    roxygenize('valley')
    R CMD INSTALL valley


