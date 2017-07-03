#!/bin/bash
echo "Files for the following runs do not exist:"
while read RUN; do
  if [ ! -f "/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED/Analysis_r0000${RUN}_s000_p000.root" ]; then
    echo $RUN
  fi
done <$1
