CellRangerSummary = FALSE
AntibodyExplore = FALSE
QCInit = FALSE
SCTxform = TRUE
ExploreBiotypes = FALSE

knitr::opts_knit$set(progress = FALSE)

### Summarize CellRanger Output Metrics
print("CellRangerSummary")
if (CellRangerSummary) {
  rmarkdown::render("0.AggregateCellRangerStats.Rmd", output_file = "../results/0.AggregateCellRangerStats.html")
}

### Explore influence of normalization and scaling on HTO and ADT outputs
print("AntibodyExplore chunk")
if (AntibodyExplore) {
  rmarkdown::render("1a.EvalMultiSeqDemux.Rmd", output_file = "../results/1a.EvalMultiSeqDemux.html")
}

### Per-run QC and initialization
print("QCInit chunk")
# List of runs that will be processed
if (QCInit) {
  runs2process = c("wk04_mpssr_wt", "wk04_mpssr_mut", 
                   "wk12_mpssr_wt", "wk12_mpssr_mut",
                   "wk36_mpssr_wt", "wk36_mpssr_mut", "wk36_novo_wt", "wk36_novo_mut")
  ADTskip = c("wk04_mpssr_mut", "wk12_mpssr_wt", "wk12_mpssr_mut", "wk36_mpssr_mut")
  
  for (myrun in runs2process) { 
    print(sprintf("Processing %s", myrun))
    calcADT = TRUE
    if (myrun %chin% ADTskip) { calcADT = FALSE }
    rmarkdown::render("/Users/lusardi/Documents/CEDAR/Projects/3.ASXL_mutant/analysis/code/1.RunQC.Rmd",
                      output_file = sprintf("../results/%s_%s.html", myrun, Sys.Date()), 
                      params = list(run2process = myrun,
                                    file2save = myrun,
                                    calcADT = calcADT))
  }
}

### Run SCTransform on all files
if (SCTxform) {
  runs2process = c("wk04_mpssr_wt", "wk04_mpssr_mut", 
                   "wk12_mpssr_wt", "wk12_mpssr_mut",
                   "wk36_mpssr_wt", "wk36_mpssr_mut", "wk36_novo_wt", "wk36_novo_mut")
  
  for (myrun in runs2process[3:8]) { 
    print(sprintf("Processing %s", myrun))
    rmarkdown::render("/Users/lusardi/Documents/CEDAR/Projects/3.ASXL_mutant/analysis/code/2a.SCTransform.Rmd",
                      output_file = sprintf("../results/2a.SCTxform_%s_%s.html", myrun, Sys.Date()), 
                      params = list(run2process = myrun,
                                    file2save = myrun,
                                    saveObj = TRUE))
  }
}
