library(dplyr)
library(tidyr)

InputPath <- "/Users/vsubasri/Medical Biophysics 2017:2018/Rotation Project 3/LFS-SNVs/"
OutputPath <- "/Users/vsubasri/Medical Biophysics 2017:2018/Rotation Project 3/SNVs_bed/"

myFiles <- list.files(path=InputPath,  pattern=".rda", recursive="TRUE")

for(fileIter in 1:length(myFiles)){ 
  input = paste(InputPath, myFiles[fileIter], sep='')
  output = paste(OutputPath, sub('.rda', '.Bed', basename(myFiles[fileIter])), sep='')
  data <- get(load(input)) 
  condensed <- data[c(1:5,38)]
  write.table(condensed, file=output, quote = FALSE, sep = "\t", row.names = FALSE)
  message(paste('File: ', output))
}
