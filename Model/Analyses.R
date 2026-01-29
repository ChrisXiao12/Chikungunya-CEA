##*****************************************************************************************************
## ANALYSES
##*****************************************************************************************************

## Base case cea
source(file = file.path(model_dir, "BaseCaseAnalysis.R"))

## One-way sensitivity analysis and tornado plots
source(file = file.path(model_dir, "OneWaySensitivityAnalysis.R"))

## Additional sensitivity analyses
source(file = file.path(model_dir, "AdditionalSensitivityAnalyses.R"))

## Probabilistic Sensitivity Analyses (PSA)
source(file = file.path(model_dir, "ProbabilisticSensitivityAnalyses.R"))

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
  