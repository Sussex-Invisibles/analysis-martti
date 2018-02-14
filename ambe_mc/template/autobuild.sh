#!/bin/bash
if [ "$#" -ne 2 ]; then
  echo "ERROR: Need exactly 2 arguments! Material, number of jobs"
  return 1
fi
export RAT_SOURCE=/home/m/mn/mn372/Software/SNOP/env_production.sh

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
export SRCDIR=/home/m/mn/mn372/Software/SNOP/work/analysis/ambe_mc/template
if [ ! -d "$SRCDIR" ]; then
  echo "ERROR: Source directory does not exist: ${SRCDIR}"
  return 3
fi
echo "Using templates from ${SRCDIR}..."

# Check for target directory
export ENDDIR=/home/m/mn/mn372/Software/SNOP/work/analysis/ambe_mc/${1}
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
  mkdir $ENDDIR/macros
  mkdir $ENDDIR/output
  mkdir $ENDDIR/results
elif [ "$(du -s $ENDDIR | awk '{print $1;}')" -ne 0 ]; then
  echo "ERROR: Target directory not empty!"
  return 5
fi

# Copy macros to target directory
echo "Copying templates to target directory..."
cp ${SRCDIR}/AmBe_TEMPLATE.geo ${ENDDIR}/AmBe_${1}.geo
cp ${SRCDIR}/AmBe_TEMPLATE.job ${ENDDIR}/AmBe_${1}.job
cp ${SRCDIR}/AmBe_TEMPLATE.mac ${ENDDIR}/AmBe_${1}.mac

# Edit macros in place
echo "Applying SED magic..."
sed -i 's|TEMPLATE|'${1}'|g' ${ENDDIR}/AmBe_${1}.geo
sed -i 's|DIRECTORY|'"${ENDDIR}"'|g' ${ENDDIR}/AmBe_${1}.job
sed -i 's|RAT_SOURCE_SCRIPT|'"${RAT_SOURCE}"'|g' ${ENDDIR}/AmBe_${1}.job
sed -i 's|TEMPLATE|'${1}'|g' ${ENDDIR}/AmBe_${1}.job
sed -i 's|DIRECTORY|'"${ENDDIR}"'|g' ${ENDDIR}/AmBe_${1}.mac
sed -i 's|TEMPLATE|'${1}'|g' ${ENDDIR}/AmBe_${1}.mac

# Copy geometry file to RAT data directory
echo "Copying geometry file..."
cp ${ENDDIR}/AmBe_${1}.geo ${RATROOT}/data/geo/calib/

# Create macros to run on cluster
echo "Generating macros..."
python ${SRCDIR}/create_macros.py ${ENDDIR} ${1} ${2}

# Submit jobs
echo "Submitting ${2} jobs..."
cd $ENDDIR
module load sge
qsub -jc mps.medium -q mps.q -t 1-${2} AmBe_${1}.job
echo "INFO: Simulating `expr 1000 \* ${2}` events for AmBeSource with ${1} as neutron absorber."
cd $SRCDIR

# Unset environment variables
unset RAT_SOURCE
unset SRCDIR
unset ENDDIR

