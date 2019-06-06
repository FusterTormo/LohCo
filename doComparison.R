tb <- read.table("logRcomp.tsv", sep = "\t", header = TRUE)
ran <- c(min(tb), max(tb))
plot(tb, xlim = ran, ylim = ran, pch = 16)
lines(x = ran, y = ran, col = "red")