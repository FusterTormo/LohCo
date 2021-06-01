#######################################################################
# Checking the combination of LOH + pathogenic variant in {gene} gene #
#######################################################################

{params}

# Load the data
{positive}
{negative}
{unknown}

# Calculate the LOH percentage reported in {gene} in positive cases
loha <- 100 * sum(pa[2:3])/sum(pa)
lohf <- 100 * sum(pf[2:3])/sum(pf)
lohn <- 100 * sum(pn[2:3])/sum(pn)
lohp <- 100 * sum(pp[2:3])/sum(pp)
lohs <- 100 * sum(ps[2:3])/sum(ps)
cat("INFO: LOH (CNN-LOH + CNL-LOH) detected in positive cases\n")
cat("\tASCAT2 reported ", round(loha,2), "% LOH in positive cases\n")
cat("\tFACETS reported ", round(lohf,2), "% LOH in positive cases\n")
cat("\tascatNGS reported ", round(lohn,2), "% LOH in positive cases\n")
cat("\tPURPLE reported ", round(lohp,2), "% LOH in positive cases\n")
cat("\tSequenza reported ", round(lohs,2), "% LOH in positive cases\n")

# Calculate the LOH percentage reported in {gene} in negative cases
loha <- 100 * sum(na[2:3])/sum(na)
lohf <- 100 * sum(nf[2:3])/sum(nf)
lohn <- 100 * sum(nn[2:3])/sum(nn)
lohp <- 100 * sum(np[2:3])/sum(np)
lohs <- 100 * sum(ns[2:3])/sum(ns)
cat("INFO: LOH detected in negative cases\n")
cat("\tASCAT2 reported ", round(loha,2), "% LOH in negative cases\n")
cat("\tFACETS reported ", round(lohf,2), "% LOH in negative cases\n")
cat("\tascatNGS reported ", round(lohn,2), "% LOH in negative cases\n")
cat("\tPURPLE reported ", round(lohp,2), "% LOH in negative cases\n")
cat("\tSequenza reported ", round(lohs,2), "% LOH in negative cases\n")

# Calculate the LOH percentage reported in {gene} in unknown cases
loha <- 100 * sum(ua[2:3])/sum(ua)
lohf <- 100 * sum(uf[2:3])/sum(uf)
lohn <- 100 * sum(un[2:3])/sum(un)
lohp <- 100 * sum(up[2:3])/sum(up)
lohs <- 100 * sum(us[2:3])/sum(us)
cat("INFO: LOH detected in unknown cases\n")
cat("\tASCAT2 reported ", round(loha,2), "% LOH in unknown cases\n")
cat("\tFACETS reported ", round(lohf,2), "% LOH in unknown cases\n")
cat("\tascatNGS reported ", round(lohn,2), "% LOH in unknown cases\n")
cat("\tPURPLE reported ", round(lohp,2), "% LOH in unknown cases\n")
cat("\tSequenza reported ", round(lohs,2), "% LOH in unknown cases\n")

# Get the LOH percentage without considering the variant classification
ma <- matrix(c(pa, na, ua), nrow = 3, byrow = TRUE)
mf <- matrix(c(pf, nf, uf), nrow = 3, byrow = TRUE)
mn <- matrix(c(pn, nn, un), nrow = 3, byrow = TRUE)
mp <- matrix(c(pp, np, up), nrow = 3, byrow = TRUE)
ms <- matrix(c(ps, ns, us), nrow = 3, byrow = TRUE)
loha <- 100 * sum(ma[,2], ma[,3]) / sum(ma)
lohf <- 100 * sum(mf[,2], mf[,3]) / sum(mf)
lohn <- 100 * sum(mn[,2], mn[,3]) / sum(mn)
lohp <- 100 * sum(mp[,2], mp[,3]) / sum(mp)
lohs <- 100 * sum(ms[,2], ms[,3]) / sum(ms)
cat("INFO: Without considering the variant classification\n")
cat("\tASCAT2 reported ", round(loha,2), "% LOH in all cases\n")
cat("\tFACETS reported ", round(lohf,2), "% LOH in all cases\n")
cat("\tascatNGS reported ", round(lohn,2), "% LOH in all cases\n")
cat("\tPURPLE reported ", round(lohp,2), "% LOH in all cases\n")
cat("\tSequenza reported ", round(lohs,2), "% LOH in all cases\n")

# Create the barplots
top <- max(ma, mf, mn, mp, ms)
colors <- c("red", "green", "yellow")
xtag <- c("Amp", "LOH", "Del", "Norm", "Not_Av")
png("ASCAT2_abs_{gene}.png", width = 730, height = 554)
barplot(ma, beside = TRUE, col = colors, names.arg = xtag, main = "ASCAT2 aberrations reported in {gene}", ylim = c(0,top))
dev.off()
png("FACETS_abs_{gene}.png", width = 730, height = 554)
barplot(mf, beside = TRUE, col = colors, names.arg = xtag, main = "FACETS aberrations reported in {gene}", ylim = c(0,top))
dev.off()
png("ascatNGS_abs_{gene}.png", width = 730, height = 554)
barplot(m, beside = TRUE, col = colors, names.arg = xtag, main = "ascatNGS aberrations reported in {gene}", ylim = c(0,top))
dev.off()
png("PURPLE_abs_{gene}.png", width = 730, height = 554)
barplot(mp, beside = TRUE, col = colors, names.arg = xtag, main = "PURPLE aberrations reported in {gene}", ylim = c(0,top))
dev.off()
png("Sequenza_abs_{gene}.png", width = 730, height = 554)
barplot(ma, beside = TRUE, col = colors, names.arg = xtag, main = "Sequenza aberrations reported in {gene}", ylim = c(0,top))
dev.off()

# Convert the absolute values to percentages
ma[1,] <- ma[1,]/sum(ma[1,])*100
ma[2,] <- ma[2,]/sum(ma[2,])*100
ma[3,] <- ma[3,]/sum(ma[3,])*100
mf[1,] <- mf[1,]/sum(mf[1,])*100
mf[2,] <- mf[2,]/sum(mf[2,])*100
mf[3,] <- mf[3,]/sum(mf[3,])*100
mn[1,] <- mn[1,]/sum(mn[1,])*100
mn[2,] <- mn[2,]/sum(mn[2,])*100
mn[3,] <- mn[3,]/sum(mn[3,])*100
mp[1,] <- mp[1,]/sum(mp[1,])*100
mp[2,] <- mp[2,]/sum(mp[2,])*100
mp[3,] <- mp[3,]/sum(mp[3,])*100
ms[1,] <- ms[1,]/sum(ms[1,])*100
ms[2,] <- ms[2,]/sum(ms[2,])*100
ms[3,] <- ms[3,]/sum(ms[3,])*100

# Create the new plots using the percentages rather than the absolute numbers
png("ASCAT2_abs_{gene}.png", width = 730, height = 554)
barplot(ma, beside = TRUE, col = colors, names.arg = xtag, main = "ASCAT2 aberrations reported in {gene}", ylim = c(0,100))
dev.off()
png("FACETS_abs_{gene}.png", width = 730, height = 554)
barplot(mf, beside = TRUE, col = colors, names.arg = xtag, main = "FACETS aberrations reported in {gene}", ylim = c(0,100))
dev.off()
png("ascatNGS_abs_{gene}.png", width = 730, height = 554)
barplot(m, beside = TRUE, col = colors, names.arg = xtag, main = "ascatNGS aberrations reported in {gene}", ylim = c(0,100))
dev.off()
png("PURPLE_abs_{gene}.png", width = 730, height = 554)
barplot(mp, beside = TRUE, col = colors, names.arg = xtag, main = "PURPLE aberrations reported in {gene}", ylim = c(0,100))
dev.off()
png("Sequenza_abs_{gene}.png", width = 730, height = 554)
barplot(ma, beside = TRUE, col = colors, names.arg = xtag, main = "Sequenza aberrations reported in {gene}", ylim = c(0,100))
dev.off()

# keys needed: gene, positive, negative, unknown, params