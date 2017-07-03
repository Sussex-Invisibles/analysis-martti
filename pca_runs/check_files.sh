#!/bin/bash
echo "Files for the following runs do not exist:"
export FPATH="/lustre/scratch/epp/neutrino/snoplus/TELLIE_PCA_RUNS_PROCESSED"
while read RUN; do
  if [ ! -f "${FPATH}/Analysis_r0000${RUN}_s000_p000.root" ] && [ ! -f "${FPATH}/Analysis_r0000${RUN}_s000_p001.root" ] && [ ! -f "${FPATH}/Analysis_r0000${RUN}_s000_p002.root" ] && [ ! -f "${FPATH}/Analysis_r0000${RUN}_s000_p003.root" ]; then
    echo $RUN
  fi
done <$1
