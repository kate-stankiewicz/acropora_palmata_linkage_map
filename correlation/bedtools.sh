for win in 1mb 25kb 40kb; do
	for chr in $(cat chroms.txt); do
		for type in palm cerv rec; do
			bedtools map -a ./${win}/ref/${chr}_filt_${win}.bed -b ./${type}/${chr}_${type}.bed -c 4 -o mean > ./${win}/map/${chr}_${win}_${type}.txt
		done
	done
done
