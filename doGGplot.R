library(ggplot2)

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
gg <- read.table("tab4ggplot_TCN.tsv", sep = "\t", header = TRUE)
gg$chr <- as.numeric(gg$chr)
ggplot(data=gg, aes(x=start, xend=end, y=cn, yend=cn, color=type)) + geom_segment() + facet_grid(.~chr) + coord_cartesian(ylim = c(0,7)) + 
  labs(title="Total copy number", x="Chromosomal_positions", y="Total_copy_number_counts")
ggsave("Total_CN.png", width = 8.27) #DIN-A4 format

gg2 <- read.table("tab4ggplot_lcn.tsv", sep = "\t", header = TRUE)
gg2$chr <- as.numeric(gg2$chr)
ggplot(data=gg2, aes(x=start, xend=end, y=cn, yend=cn, color=type)) + geom_segment() + facet_grid(.~chr) + coord_cartesian(ylim = c(0,7)) + 
  labs(title="Minor copy number", x="Chromosomal_positions", y="Minor_copy_number_counts")
ggsave("Minor_CN.png", width = 8.27)