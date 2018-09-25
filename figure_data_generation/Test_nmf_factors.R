#
#
#
#
#
#

rm(list = ls())

library(NMF)
setwd("E:/Ori_F/ZWen/TestNMF")

# Create data matrix
padata <- read.csv("fusion_matrix.csv", header = TRUE)

res <- nmf(padata, 2:16, .opt='vp2', nrun=50)

rt <- randomize(padata)
res_rt <- nmf(rt, 2:16, .opt='vp2', nrun=50)

pdf("final_run50.pdf")
plot(res, res_rt)
dev.off()