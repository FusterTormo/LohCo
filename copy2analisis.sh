#! /bin/bash

counter=0

cd ~/Desktop/compareAscatFacets

if [ ! -d analisis ];then
	mkdir analisis
else
	rm analisis/*
fi

for i in */*_tab4ggplot_lcn.tsv
do
	cp $i ~/Desktop/compareAscatFacets/analisis
	temp=`basename $i`
	mv analisis/${temp} analisis/${counter}_${temp}
	counter=$((counter+1))
done



echo "Copied $counter files"
