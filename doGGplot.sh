#! /bin/bash

counter=0

cd ~/Desktop/compareAscatFacets

for i in */*_tab4ggplot_lcn.tsv
do
	d=`dirname $i`
	echo "INFO: Plotting $i"
	Rscript ~/Desktop/compareAscatFacets/LohCo/doGGplot.R $i
	echo "INFO: Moving plot to $d"
	mv Minor_CN.png $d
done
