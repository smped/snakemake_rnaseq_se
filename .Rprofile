## Ensure that the conda installations are used, if running R from a conda environment
local({
  conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
  if (!is.na(conda) && startsWith(R.home(), conda)) {
    .libPaths(R.home("library"))
  }
})

## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}
