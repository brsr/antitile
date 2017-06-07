#!/bin/bash 
#IFS='/' tokens=($1)
#coda="${tokens[-1]}"
#for a in {1..9}; do	for ((b=0; b<=$a; b++)); do
	#geodesic -c "$a","$b" "$1" -o "$a-$b-anti-n-0-0-$coda"
	#geodesic -c "$a","$b" "$1" | minmax -n 100000 -o "$a-$b-anti-m-0-0-$coda"	
	#geodesic -c "$a","$b" "$1" | repel -n 100000 -o "$a-$b-anti-r-0-0-$coda"
	#minmax -n 100000 "0$a-0$b-flat-n-0-0-cube.off" -o "$a-$b-anti-m-0-0-cube.off"
	#repel -n 100000 "0$a-0$b-flat-n-0-0-cube.off" -o  "$a-$b-anti-r-0-0-cube.off"
#done; done

#for filename in  *-nslerp-n-1-1-*; do
	#outfile=${filename/nslerp-n/anti-c}
	#off_util -D vo2 $filename | canonical -o $outfile
#done
