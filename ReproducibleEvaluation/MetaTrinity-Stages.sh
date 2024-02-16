#!/bin/sh


echo "\n=======================================================  RUNNING STAGE 1: Containment Search  ======================================================\n"

/usr/bin/time -v --output Results/ContainmentSearch/ContainmentSearch_time.txt python3 ../MetaTrinity/Scripts/containment_search.py ../Data/ReadSets/test.fastq ../Data/RefData/test --mmi_dir ../Data/RefData/test --temp_dir Results/ContainmentSearch --translation ../Data/RefData/test/translate/translate_sorted.csv --keep_temp_files

echo "\n=======================================================  DONE STAGE 1 ======================================================"



echo "\n=======================================================  RUNNING STAGE 2: Read Mapping  ======================================================\n"

./ReadMapping.sh -s 15 -e 15 -i 1 -n 3 -r ../Data/ReadSets/test.fastq -d Results/ContainmentSearch -t Results/ReadMapping/Timing -o Results/ReadMapping/SAM-files

echo "\n=======================================================  DONE STAGE 2 ======================================================"




echo "\n=======================================================  RUNNING STAGE 3: Taxonomic Profiling  ======================================================\n"

./TaxonomicProfiling.sh -r 15 -n 3 -s Results/ReadMapping/SAM-files -o Results/TaxonomicProfiling

echo "\n=======================================================  DONE STAGE 3 ======================================================"

