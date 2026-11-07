# Paths on Fir. Source this first in every session. For another Alliance cluster, this is the file to edit.
# Everything heavy lives on scratch. The repo (on /project) only receives the profiler outputs at the end.
export ILAFORES=$(realpath ~/scratch)
export REPO=/project/6047190/ronj2303/MethodsComparison # git repo; data/Horse gets the profiler outputs at the end
export MC=$ILAFORES/mc_work                               # working copy: data/ on scratch, scripts -> $REPO/scripts
export ILL_PIPELINES=$ILAFORES/ILL_pipelines              # branch kneaddata_gzip
export DBS=$ILAFORES/ref_dbs                              # databases, host indexes included
export RAW=$ILAFORES/Horse/raw                            # raw reads (~550 GB). Scratch purges after ~60 days
