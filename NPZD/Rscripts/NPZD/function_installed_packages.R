 # functions to load packages that are needed for NPZ model (Steven Pint)

install_needed_packages <- function(){ 
not_installed <- setdiff(c("ggplot2", 
                           "xts", 
                           "reshape2", 
                           "lubridate", 
                           "stats", 
                           "plyr", 
                           "dplyr", 
                           "RColorBrewer", 
                           "viridis", 
                           "ggpubr", 
                           "parallel", 
                           "doParallel", 
                           "foreach", 
                           "knitr",
                           "solrad",
                           "remotes",
                           "mgcv")
                         , rownames(installed.packages()))


not_installed_github <- setdiff("lwdataexplorer", rownames(installed.packages()))

 if (length(not_installed) != 0 | length(not_installed_github) != 0) {
  
   if (length(not_installed) != 0) {
     install.packages(not_installed)
     }
   if (length(not_installed_github) != 0) {
     remotes::install_github("lifewatch/lwdataexplorer")
   }  
   } else {
   print("All packages that are needed are installed")
 }
}

load_all_packages <- function(){
    all_packages <- c("ggplot2", 
                           "xts", 
                           "reshape2", 
                           "lubridate", 
                           "stats", 
                           "plyr", 
                           "dplyr", 
                           "RColorBrewer", 
                           "viridis", 
                           "ggpubr", 
                           "parallel", 
                           "doParallel", 
                           "foreach", 
                           "knitr",
                           "solrad",
                           "remotes",
                           "mgcv",
                     "lwdataexplorer")
lapply(all_packages, require, character.only = TRUE)
}