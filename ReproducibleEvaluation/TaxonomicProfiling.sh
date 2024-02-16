#!/bin/sh

while getopts r:n:s:o: flag
do
    case "${flag}" in
        r) r=${OPTARG};;
        n) n=${OPTARG};;
		    s) sam=${OPTARG};;
		    o) out=${OPTARG};;
    esac
done

echo "\n================================== INPUT =================================="
echo "Edit distance threshold (-r): $r";
echo "Number of mapping locations (-n): $n"
echo "SAM Files: $sam";
echo "\n================================== OUTPUT  =================================="
echo "Output directory: $out";
echo "\n================================== ALGORITHMS =================================="
echo "adjacency-filter base-counting edlib grim_original grim_original_tweak hd magnet qgram shd shouji sneakysnake \n\n"

DIR="../MetaFast"


for filter in grim_original_tweak magnet hd adjacency-filter base-counting qgram shd shouji sneakysnake grim_original edlib
  do
    echo "Profiling for "$filter" "
    /usr/bin/time -v --output "$out"/Timing/profiling-time__"$filter"_n"$n"_r"$r".txt python3 $DIR/Scripts/read_mapping.py $sam/"$filter"__n"$n"_r"$r".sam /home/arvidg/MetaFast/Data/RefData/test/ --input_type sam --verbose --output "$out"/Profiles/"$filter"__n"$n"_r"$r".tsv
    echo "================================================ \n\n"
  done
