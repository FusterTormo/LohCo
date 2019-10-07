#########################################################################
#Creates a ggplot from a table file passed as parameter.                #
#The structure of the table should be:                                  #
# chr start end CN  type                                                #
#Where  chr is the chromosome                                           #
#       start is the starting position of the segment                   #
#       end is the end position of the segment                          #
#       CN is the copy number value of the segment                      #
#       type is the type of program that has inferred this segment      #
#########################################################################
arg <- commandArgs(trailingOnly = TRUE)

if (length(arg) >= 2) {
  fil <- arg[1]
  output <- arg[2]
} else if (length(arg) == 1) {
  fil <- arg[1]
  output <- "CN.png"
} else {
  stop("ERROR: Not enough arguments passed to R script.\nUsage:\n\tRscript doGGplot.R input_file.tsv [output_file.png]")
}

library(ggplot2)

gg <- read.table(fil, sep = "\t", header = TRUE)
gg$chr <- as.numeric(gg$chr)
ggplot(data=gg, aes(x=start, xend=end, y=cn, yend=cn, color=type)) + geom_segment() + facet_grid(.~chr) + coord_cartesian(ylim = c(0,4)) + 
  labs(title="Copy number", x="Chromosomal_positions", y="Copy_number_counts")
ggsave(output, width = 8.27) #DIN-A4 format