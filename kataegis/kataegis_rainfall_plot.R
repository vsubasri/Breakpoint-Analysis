library(karyoploteR)
library(regioneR)

setwd("~/research/lfs somatic analysis/LFS-SNVs")

preprocess <- function(data){
  load(data)
  data <- data.filtered[c(1:5,38)]
  colnames(data) <- c("chr", "start", "end", "ref", "alt", "mut.type")
  sm.gr <- toGRanges(data[,c("chr", "start", "end", "mut.type", "ref", "alt")])
  seqlevelsStyle(sm.gr) <- "UCSC"
  return(sm.gr)
} 

  file <-'2821B_annotated_filtered_clipped.rda'
  title <- '2821B Kataegis'
  sm.gr <- preprocess(file)
  variant.colors <- getVariantsColors(sm.gr$ref, sm.gr$alt)
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 0
  pp$bottommargin <- 20
  kp <- plotKaryotype(plot.type=4, ideogram.plotter = kpAddCytobandsAsLine,labels.plotter = NULL, plot.params = pp)
  kpAddChromosomeNames(kp, srt=45)
  kpAddMainTitle(kp, main=title, cex=1.2)
  kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0, r1=0.65)
  kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.65)
  kpAddLabels(kp, labels = 'Distance between mutations (log10)', srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.7)
  kpPlotDensity(kp, data = sm.gr, r0=0.72, r1=1)
  kpAddLabels(kp, labels = 'Density', srt=90, pos=1, label.margin = 0.045, r0=0.71, r1=1)
  legend("topright", 
       c("C>A", "C>T", "C>G", "T>C", "T>G", "T>A"),
       col = c("#4c64ae", "#e40611", "#000000", "#fbe800", "#6eb529","#bf4a96"),
       pch = c(20),
       bty = "n",
       pt.cex = 0.8,
       cex = 0.4,
       text.col = "black",
       horiz = F,
       xpd =TRUE,
       inset=c(-0.12,0.5))

##per chromosome 
chromosome <- function(file, chr, title) {
  sm.gr <- preprocess(file)
  pp <- getDefaultPlotParams(plot.type = 4)
  variant.colors <- getVariantsColors(sm.gr$ref, sm.gr$alt)
  kp <- plotKaryotype(plot.type=4, plot.params = pp, chromosomes=chr)
  kpAddMainTitle(kp, main=title, cex=1.2)
  kpAddBaseNumbers(kp)
  kpPlotDensity(kp, data = sm.gr, window.size = 10e5, r0=0.62, r1=0.8)
  kpAddLabels(kp, labels = 'Density', srt=90, pos=1, label.margin = 0.048, r0=0.62, r1=0.8)
  kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0, r1=0.6)
  kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.6)
  kpAddLabels(kp, labels = 'Distance between mutations (log10)', srt=90, pos=1, label.margin = 0.048, r0=0, r1=0.6)
  legend("topright", 
       c("C>A", "C>T", "C>G", "T>C", "T>G", "T>A"),
       col = c("#4c64ae", "#e40611", "#000000", "#fbe800", "#6eb529","#bf4a96"),
       pch = c(20),
       bty = "n",
       pt.cex = 0.8,
       cex = 0.4,
       text.col = "black",
       horiz = F,
       xpd =TRUE,
       inset=c(-0.12,0.5))
}

for (i in 17:length(files_snv)) {
  rainfall(files_snv[i], basename(files_snv[i]))
}
chromosome('4856_3_annotated_filtered_clipped.rda', 'chr11', 'Somatic Mutations - 4333 Chr11') 
