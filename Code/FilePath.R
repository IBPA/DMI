PATH = dirname(this.path::this.dir())
DATA_PATH = paste0(PATH,"\\Data")
CODE_PATH = paste0(PATH,"\\Code")
setwd(PATH)

PROCESSED_DATA_PATH = paste0(PATH,"\\ProcessedData")
FIGURE_PATH=paste0(PATH,"\\Figures")

dir.create(PROCESSED_DATA_PATH, showWarnings = FALSE)
dir.create(FIGURE_PATH, showWarnings = FALSE)


