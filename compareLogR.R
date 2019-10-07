sysarg <- commandArgs(trailingOnly = TRUE)
fil <- sysarg[1][1]
tb <- read.table(fil, sep = "\t", header = TRUE)

ran <- c(min(tb), max(tb))
png("logRcorrelation.png", height = 700, width = 700)
plot(tb, xlim = ran, ylim = ran, pch = 16)
lines(x = ran, y = ran, col = "gray")
dev.off()