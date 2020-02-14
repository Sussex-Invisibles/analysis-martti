#!/bin/bash
echo "Files for the following runs do not exist:"
export FPATH="/its/home/mr514/DEC_2019/calib"
#export FPATH="/lustre/scratch/epp/neutrino/snoplus/TELLIE/TELLIE_PCA_RUNS"
#export FPATH="/lustre/scratch/epp/neutrino/snoplus/TELLIE/TELLIE_STABILITY_RUNS_PROCESSED"
#export FPATH="$HOME/Software/SNOP/work/data"
#export FPATH="downloaded"
while read RUN; do
  for f in ${FPATH}/*${RUN}*; do
    if [ ! -e $f ]; then
      echo $RUN
    fi
    break
  done
done <$1
