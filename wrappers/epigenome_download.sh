#!/bin/bash

declare -a ipsc_epis=(
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E018_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E019_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E020_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E021_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E022_15_coreMarks_hg38lift_mnemonics.bed.gz"
  )
  
declare -a heart_epis=(
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E083_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E104_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E095_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E105_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E065_15_coreMarks_hg38lift_mnemonics.bed.gz"
  )

declare -a muscle_epis=(
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E078_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E076_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E103_15_coreMarks_hg38lift_mnemonics.bed.gz"
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E111_15_coreMarks_hg38lift_mnemonics.bed.gz"
  )

for f in ${ipsc_epis[@]}; do
  fname=`grep -oP '(?<=final/).*?(?=_)' <<< "$f"`
  wget -O "../data/epigenomes/iPSC/"$fname"_15state.bed.gz" $f 
  gzip -d "../data/epigenomes/iPSC/"$fname"_15state.bed.gz"
done

for f in ${heart_epis[@]}; do
  fname=`grep -oP '(?<=final/).*?(?=_)' <<< "$f"`
  wget -O "../data/epigenomes/heart/"$fname"_15state.bed.gz" $f 
  gzip -d "../data/epigenomes/heart/"$fname"_15state.bed.gz"
done

for f in ${muscle_epis[@]}; do
  fname=`grep -oP '(?<=final/).*?(?=_)' <<< "$f"`
  wget -O "../data/epigenomes/muscle/"$fname"_15state.bed.gz" $f 
  gzip -d "../data/epigenomes/muscle/"$fname"_15state.bed.gz"
done
