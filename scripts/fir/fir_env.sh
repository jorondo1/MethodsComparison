# Paths on Fir. Source this first in every session. For another Alliance cluster, this is the file to edit.
export ILAFORES=$(realpath ~/scratch) # everything runs on scratch
export ILL_PIPELINES=$ILAFORES/programs/ILL_pipelines
export FAST=$ILAFORES/fast # MetaPhlAn databases and host indexes, so that existing paths still work
export MC=$ILAFORES/MethodsComparison
export RAW=/project/def-ilafores/ronj2303/Horse/raw         # raw reads: project is not purged
export RESULTS=/project/def-ilafores/ronj2303/Horse/results # final outputs, copied from scratch

mkdir -p $RESULTS