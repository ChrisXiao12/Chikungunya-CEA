##*****************************************************************************************************
## ANALYSES
##*****************************************************************************************************

## Base case cea
if(isTRUE(run_base_case)){
  source(file = file.path(model_dir, "BaseCaseAnalysis.R"))
}

## One-way sensitivity analysis and tornado plots
if(isTRUE(run_owsa)){
  source(file = file.path(model_dir, "OneWaySensitivityAnalysis.R"))
}

## Additional sensitivity analyses
if(isTRUE(run_addi_sensi_analysis)){
  source(file = file.path(model_dir, "AdditionalSensitivityAnalyses.R"))
}

## Probabilistic Sensitivity Analyses (PSA)
if(isTRUE(run_psa)){
  source(file = file.path(model_dir, "ProbabilisticSensitivityAnalyses.R"))
}

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
  