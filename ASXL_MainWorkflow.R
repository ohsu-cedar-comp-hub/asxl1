CellRangerSummary = FALSE
AntibodyExplore = FALSE
QCInit = TRUE
SCTxform = FALSE
ExploreBiotypes = FALSE
local = FALSE
message(sprintf("local = %s", local))
if (local) {
  dir = list(code = "/Users/lusardi/Documents/CEDAR/Projects/3.ASXL_mutant/analysis/code",
                   results = "/Users/lusardi/Documents/CEDAR/Projects/3.ASXL_mutant/analysis/results")
} else {
  dir  = list(code = "/home/groups/CEDAR/lusardi/asxl/asxl1",
                   results = "/home/groups/CEDAR/braun/seurat_obj/test_rmd")
}

knitr::opts_knit$set(progress = FALSE)
library(data.table)
library(rmarkdown)
library(Seurat)

processDate = Sys.Date()

### Summarize CellRanger Output Metrics
print("CellRangerSummary")
if (CellRangerSummary) {
  rmarkdown::render("0.AggregateCellRangerStats.Rmd",
                    output_file = sprintf("%s/0.AggregateCellRangerStats_%s.html", dir$results, processDate))
}

### Explore influence of normalization and scaling on HTO and ADT outputs
print("AntibodyExplore chunk")
if (AntibodyExplore) {
  rmarkdown::render("1a.EvalMultiSeqDemux.Rmd",
                    output_file = sprintf("%s/1a.EvalMultiSeqCDemux_%s.html", dir$results, processDate))
}

### Per-run QC and initialization
print("QCInit chunk")
if (QCInit) {
  # List of runs that will be processed
  runs2process = c("wk04_mpssr_wt", "wk04_mpssr_mut", 
                   "wk12_mpssr_wt", "wk12_mpssr_mut",
                   "wk36_mpssr_wt", "wk36_mpssr_mut", "wk36_novo_wt", "wk36_novo_mut")
  ADTskip = c("wk04_mpssr_mut", "wk12_mpssr_wt", "wk12_mpssr_mut", "wk36_mpssr_mut")
  sourceDir = "/home/groups/CEDAR/braun/cellranger"

  # Set up destination directory
  baseDir = "/home/groups/CEDAR/braun/seurat_obj"
  destDir = sprintf("%s/1.RunQC_%s", baseDir, processDate)
  # If destDir does not exist, create it
  if (!dir.exists(destDir)) {
      dir.create(destDir)
    if (!dir.exists(destDir)) {
        message(sprintf("ERROR: could not create directory %s", destDir))
        knitr::knit_exit()
    }
  }

  for (myrun in runs2process) { 
    print(sprintf("Processing input run %s", myrun))
    calcADT = TRUE
    if (myrun %chin% ADTskip) { calcADT = FALSE }
    print(sprintf("Processing script %s", paste(dir$code, "1.RunQC.Rmd", sep = "/")))
    print(dir$code)
    rmarkdown::render(paste(dir$code, "1.RunQC.Rmd", sep = "/"),
                      output_file = sprintf("%s/1.%s_%s.html", destDir, myrun, processDate),
                      params = list(run2process = myrun,
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
                      output_file = sprintf("../results/2a.SCTxform_%s_%s.html", myrun, processDate), 
                      params = list(run2process = myrun,
                                    file2save = myrun,
                                    saveObj = TRUE))
  }
}
