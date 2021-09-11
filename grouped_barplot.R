library(ggplot2)
library(reshape)

dat <- read.csv('snvs_breakpoint_1kb_stats.tsv', sep='\t', col.names=c("sample", "count", "mean", "std", "close", "far", "TA","TG", "TC", "CT","CA","CG"))
dat[is.na(dat)] <- 0
dat <- dat[c(1,7:12)]
data.m <- melt(dat, id.vars='sample')
#blue red black yellow green purple
df = data.frame(muts=c("CA", "CT", "CG", "TC", "TG", "TA"), color=c("#4c64ae", "#e40611", "#000000", "#fbe800", "#6eb529","#bf4a96"))
data.m["color"] <- df$color[match(data.m$variable,df$mut)]
ggplot(data.m, aes(sample, value, fill=variable)) + 
  geom_bar(position = "dodge", stat="identity") +
  theme_bw() +
  theme(legend.key.height=unit(0.5,"cm")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))+
  labs(x="Sample",y="Number of SNVs within 1kb of Breakpoints",fill="Mutation Type") + 
  scale_fill_manual(name = "",
                    labels=c("TA", "TG", "TC", "CT", "CA", "CG"),
                    values=c("#bf4a96", "#6eb529", "#fbe800", "#e40611", "#4c64ae","#000000"))
  