#!/bin/sh

while getopts s:e:i:n:r:d:t:o: flag
do
    case "${flag}" in
        s) r_start=${OPTARG};;
        e) r_end=${OPTARG};;
        i) r_increment=${OPTARG};;
        n) n=${OPTARG};;
		    r) reads=${OPTARG};;
		    d) subset_DB=${OPTARG};;
        t) timing=${OPTARG};;
        o) out=${OPTARG};;
    esac
done
echo "\n================================== INPUT =================================="
echo "Start value for -r: $r_start";
echo "End value for -r: $r_end";
echo "Increment for -r: $r_increment";
echo "Number of mapping locations per read: $n"
echo "Reads: $reads";
echo "Subset DB: $subset_DB";
echo "\n================================== OUTPUT  =================================="
echo "Results/ReadMapping";
echo "\n================================== ALGORITHMS =================================="
echo "adjacency-filter base-counting edlib grim_original grim_original_tweak hd magnet qgram shd shouji sneakysnake"

DIR="../MetaFast/ReadMapping"


for filter in hd adjacency-filter base-counting magnet qgram shd shouji sneakysnake edlib grim_original grim_original_tweak
  do
    for r in $(seq $r_start $r_increment $r_end)
      do
        echo "Read Mapping with "$filter" "
        /usr/bin/time -v --output $timing/time__"$filter"_n"$n"_r"$r".txt $DIR/rm -ax sr -t 1 -2 -n $n -r $r --filter="$filter" --sam-hit-only --secondary=yes $subset_DB/subset_db.fna $reads > $out/"$filter"__n"$n"_r"$r".sam
        echo "================================================ \n"
      done
  done