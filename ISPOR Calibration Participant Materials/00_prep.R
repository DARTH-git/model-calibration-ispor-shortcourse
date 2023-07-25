### Run all code in this script to download the desired R packages 
### from either CRAN or GitHub


# first download and use this package to conveniently install other packages
if (!require('pacman')) install.packages('pacman'); library(pacman) 

# load (install if required) packages from CRAN
p_load("devtools", "lhs", "plotrix", "psych", "DEoptim", "matrixStats", "scatterplot3d")

# Install IMIS
devtools::install_version("IMIS", version = "0.1", repos = "http://cran.us.r-project.org")



          
  
          
