<p align="center">
  <img width="424" height="200" src="https://github.com/tgac-vumc/Statescope/blob/master/img/Logo_Statescope.png">
</p>

# StatescopeR
StatescopeR is an R wrapper around Statescope, a computational framework designed to discover cell states from cell type-specific gene expression profiles inferred from bulk RNA profiles.
See https://github.com/tgac-vumc/Statescope for the original Python version.
<p align="center">
  <img width="75%" height="75%" src="https://github.com/tgac-vumc/Statescope/blob/master/img/Statescope_Overview.png">
</p>


## Installation
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("StatescopeR")
