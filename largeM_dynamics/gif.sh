#!/bin/bash
#########################################################################
#File Name: gif.sh
# Author:  Tiangang Zhou
# mail: tg_zhou@pku.edu.cn
#########################################################################

# for dir in $(find . -type d | grep '/'); do
	# cd $dir
	# find . -type f -name 'ρω*_thermal.pdf' -print0 |
	# while IFS= read -r -d '' file
	# 	do convert -verbose -density 500 -resize 800 "${file}" "${file%.*}.png"
	# done
	echo "Start generate gif in " 
	cp ρω0_equilbrium.png ρω0_thermal.png
	pwd
	convert -delay 100 -dispose previous $(ls -v ρω*_thermal.png) wormhole_eq.gif 
	echo "end"
	cd .. 
# done
echo "gif generate successfully"


	# find . -type f -name '*.pdf' -print0 |
	# while IFS= read -r -d '' file
	# 	do convert -verbose -density 500 -resize 800 "${file}" "${file%.*}.png"
	# done

