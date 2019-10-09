CellRangerSummary = FALSE
AntibodyExplore = FALSE
QCInit = TRUE
SCTxform = FALSE
ExploreBiotypes = FALSE
local = FALSE
message(sprintf("local = %s", local))
if (local) {
  directory = list(code = "/Users/lusardi/Documents/CEDAR/Projects/3.ASXL_mutant/analysis/code",
                   results = "/Users/lusardi/Documents/CEDAR/Projects/3.ASXL_mutant/analysis/results")
} else {
  directory = list(code = "/home/groups/CEDAR/lusardi/asxl/asxl1",
                   results = "/home/groups/CEDAR/braun/seurat_obj/test_rmd")
}

knitr::opts_knit$set(progress = FALSE)
library(data.table)
library(rmarkdown)
library(Seurat)

### Summarize CellRanger Output Metrics
print("CellRangerSummary")
if (CellRangerSummary) {
  rmarkdown::render("0.AggregateCellRangerStats.Rmd",
                    output_file = sprintf("%s/0.AggregateCellRangerStats_%s.html", directory$results, Sys.Date()))
}

### Explore influence of normalization and scaling on HTO and ADT outputs
print("AntibodyExplore chunk")
if (AntibodyExplore) {
  rmarkdown::render("1a.EvalMultiSeqDemux.Rmd",
                    output_file = sprintf("%s/1a.EvalMultiSeqCDemux_%s.html", directory$results, Sys.Date()))
}

### Per-run QC and initialization
print("QCInit chunk")
# List of runs that will be processed
if (QCInit) {
  runs2process = c("wk04_mpssr_wt", "wk04_mpssr_mut", 
                   "wk12_mpssr_wt", "wk12_mpssr_mut",
                   "wk36_mpssr_wt", "wk36_mpssr_mut", "wk36_novo_wt", "wk36_novo_mut")
  ADTskip = c("wk04_mpssr_mut", "wk12_mpssr_wt", "wk12_mpssr_mut", "wk36_mpssr_mut")
  baseDir = "/home/groups/CEDAR/braun"
  sourceDir = "cellranger"
  destDir = paste("1.SeuratObj_TEST", Sys.Date(), sep = "_")
  
  for (myrun in runs2process[1]) { 
    print(sprintf("Processing %s", myrun))
    calcADT = TRUE
    if (myrun %chin% ADTskip) { calcADT = FALSE }
    rmarkdown::render(paste(directory$code, "1.RunQC.Rmd", sep = "/"),
                      output_file = sprintf("%s/1.%s_%s.html", directory$results, myrun, Sys.Date()),
                      params = list(run2process = myrun,
                                    baseDir = baseDir,
                                    sourceDir = sourceDir,
                                    destDir = destDir,
                                    file2save = myrun,
                                    calcADT = calcADT,
                                    local = local))
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
