#!/bin/bash
if [ "$#" -ne 2 ]; then
  echo "ERROR: Need exactly 2 arguments! Decay mode, number of jobs"
  return 1
fi
export RAT_SOURCE=/home/m/mn/mn372/Software/SNOP/env_dinucleon.sh

# Check for RAT installation
if [ -z "$RATROOT" ]; then
  echo "INFO: Need to source RAT installation."
  source ${RAT_SOURCE}
fi
if [ ! -d "$RATROOT" ]; then
  echo "ERROR: RAT directory does not exist: ${RATROOT}"
  return 2
fi
echo "Found RAT installation: ${RATROOT}"

# Check for source directory
export SRCDIR=/home/m/mn/mn372/Software/SNOP/work/analysis/dinucleon/template
if [ ! -d "$SRCDIR" ]; then
  echo "ERROR: Source directory does not exist: ${SRCDIR}"
  return 3
fi
echo "Using templates from ${SRCDIR}..."

# Check for target directory
export ENDDIR=/home/m/mn/mn372/Software/SNOP/work/analysis/dinucleon/${1}
if [ ! -d "$ENDDIR" ]; then
  mkdir $ENDDIR
fi
if [ ! -d "$ENDDIR" ]; then
  echo "ERROR: Target directory does not exist: ${ENDDIR}"
  return 4
fi
echo "Target directory is ${ENDDIR}..."

# Build subfolders in target directory
if [ -z "$(ls -A ${ENDDIR})" ]; then
  mkdir $ENDDIR/logs
  mkdir $ENDDIR/logs/bash
  mkdir $ENDDIR/logs/rat_logs
  mkdir $ENDDIR/output
  mkdir $ENDDIR/results
elif [ "$(du -s $ENDDIR | awk '{print $1;}')" -ne 0 ]; then
  echo "ERROR: Target directory not empty!"
  return 5
fi

# Copy files to target directory
echo "Copying templates to target directory..."
cp ${SRCDIR}/Dinucleon_TEMPLATE.job ${ENDDIR}/Dinucleon_${1}.job

# Edit macros in place
echo "Applying SED magic..."
sed -i 's|DIRECTORY|'"${ENDDIR}"'|g' ${ENDDIR}/Dinucleon_${1}.job
sed -i 's|RAT_SOURCE_SCRIPT|'"${RAT_SOURCE}"'|g' ${ENDDIR}/Dinucleon_${1}.job
sed -i 's|TEMPLATE|'${1}'|g' ${ENDDIR}/Dinucleon_${1}.job

# Submit jobs
echo "Submitting ${2} jobs..."
cd $ENDDIR
module load sge
qsub -jc mps.medium -q mps.q -l h_rt=12:00:00 -t 1-${2} Dinucleon_${1}.job

echo "INFO: Simulating `expr 10000 \* ${2}` events for Dinucleon Decay (${1} mode)."
cd $SRCDIR

# Unset environment variables
unset RAT_SOURCE
unset SRCDIR
unset ENDDIR

