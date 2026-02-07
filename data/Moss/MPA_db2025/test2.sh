export EXE_PATH=$(dirname "$0")

echo "EXE_PATH=$EXE_PATH"

echo "copy"
cp ${EXE_PATH}/../containers/humann.3.6.sif $SLURM_TMPDIR/
echo "copy done"
