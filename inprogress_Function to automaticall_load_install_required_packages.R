loadLibrary <- function(lib){
  # Require returns TRUE invisibly if it was able to load package
  # If package was not able to be loaded then first check if it can be loaded with Bioconductor
  # If not, attempt to reinstall with Bioconductor and then install.packages 
  # Load package after (re)installing
  for(i in lib){
    if (!require("BiocManager")){
      install.packages("BiocManager")
      if(!require(i, character.only = TRUE)){
        BiocManager::install(i)
        if(!require(i, character.only = TRUE)){
          install.packages(i, dependencies = TRUE)
          require(i, character.only = TRUE)}}}}}