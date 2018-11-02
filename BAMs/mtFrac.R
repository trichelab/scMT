library(Rsamtools) 

mtFrac <- function(bam) {
  idx <- idxstatsBam(bam)
  subset(idx, seqnames=="chrM")$mapped / sum(idx$mapped)
}

scBams <- grep("merged", invert=TRUE, list.files(patt="bam$"), value=TRUE)

scMTreads <- sapply(scBams, mtFrac)
saveRDS(scMTreads, file="scMTreads.rds") 
summary(scMTreads)

#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0241299 0.2695761 0.3853691 0.3985105 0.5259580 0.7664648 

# plot 
library(ggplot2) 
library(ggforce)
library(ggthemes) 
mtFracs <- data.frame(mtFrac=scMTreads,
                      subject=elts(names(scMTreads)),
                      celltype=elts(names(scMTreads), "_", 2))

p <- ggplot(mtFracs, aes(x=subject, color=celltype, y=mtFrac))
p + geom_sina() + 
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
  ylab("Reads mapping to chrM") +
  theme_tufte(base_size=18) + 
  ggtitle("Mitochondrial reads in single-cell ATACseq") + 
  NULL
ggsave("mtFrac.png", type="cairo")
