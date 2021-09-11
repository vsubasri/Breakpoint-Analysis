library(dplyr)
library(tidyr)

InputPath <- "/Users/vsubasri/research/filtered_tab_files/"
OutputPath <- "/Users/vsubasri/research/filtered_tab_files_bed/"

myFiles <- list.files(path=InputPath,  pattern=".tab", recursive="TRUE")

for(fileIter in 1:length(myFiles)){ 
  input = paste(InputPath, myFiles[fileIter], sep='')
  output = paste(OutputPath, sub('.tab', '.Bed', basename(myFiles[fileIter])), sep='')
  data <- read.delim(input)
  condensed <- data[c(1:4)]
  write.table(condensed, file=output, quote = FALSE, sep = "\t", row.names = FALSE)
  message(paste('File: ', output))
}
