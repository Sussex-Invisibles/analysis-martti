#!/bin/bash
module load sge

# Check input arguments
if [ "$#" -ne 1 ]; then
  echo "ERROR: Need exactly 1 argument (Material)"
  return 1
fi

# Check for RAT installation
export RAT_SOURCE=/home/m/mn/mn372/Software/SNOP/env_production.sh
if [ -z "$RATROOT" ]; then
  echo "INFO: Need to source RAT installation."
  source ${RAT_SOURCE}
fi
if [ ! -d "$RATROOT" ]; then
  echo "ERROR: RAT directory does not exist: ${RATROOT}"
  return 2
fi
echo "Found RAT installation: ${RATROOT}"

# Check for target directory
export DIR=/home/m/mn/mn372/Software/SNOP/work/analysis/ambe_mc/${1}
if [ ! -d "$DIR" ]; then
  mkdir $DIR
fi
if [ ! -d "$DIR" ]; then
  echo "ERROR: Target directory does not exist: ${DIR}"
  return 4
fi
echo "Target directory is ${DIR}..."

# Loop over bash logs
count=0
for bashlog in ${DIR}/logs/bash/*; do

  # Skip jobs that completed successfully
  last=`tail ${bashlog} -n 1`
  if [[ ${last} == "Finished job script"* ]]; then
    continue
  fi
  
  # Reset/increase counters
  (( count++ ))
  found=0
  
  # Find bash log of failed job
  index=`echo ${bashlog} | sed 's/.*\.//'`
  if [ -f ${bashlog} ] ; then 
    (( found++ ))
  else
    echo "WARNING: Could not find BASH log: ${bashlog}"
  fi
  
  # Find RAT log of failed job
  line=`grep "Hostname: node" ${bashlog}`
  node=`echo $line | cut -d " " -f 2 | sed 's/[^0-9]*//g'`
  pid=`echo $line | cut -d " " -f 4`
  ratlog=${DIR}/logs/rat_logs/rat.node${node}.${pid}.log
  if [ -f ${ratlog} ] ; then
    (( found++ ))
  else
    echo "WARNING: Could not find RAT log: ${ratlog}"
  fi
  
  # Find ROOT output file of failed job
  rootfile=${DIR}/output/AmBe_${1}_${index}.root
  if [ -f ${rootfile} ] ; then
    (( found++ ))
  else
    echo "WARNING: Could not find ROOT file: ${rootfile}"
  fi
  
  # Find macro that submitted failed job
  macro=${DIR}/macros/AmBe_${1}_${index}.mac
  if [ -f ${macro} ] ; then
    (( found++ ))
  else
    echo "WARNING: Could not find RAT macro: ${macro}"
  fi
  
  # Delete files ONLY if all 4 files were found
  if [ $found -ne 4 ]; then
    echo "WARNING: Missing output files ($found) for job #${index}. Will not delete anything."
  else
    echo "Deleting files for job #${index}."
    rm ${bashlog}
    rm ${ratlog}
    rm ${rootfile}
  fi
  
  # Resubmit jobs
  echo "qsub -jc mps.medium -q mps.q -l h_rt=12:00:00 -t ${index} ${DIR}/AmBe_${1}.job"
  qsub -jc mps.medium -q mps.q -l h_rt=12:00:00 -t ${index} ${DIR}/AmBe_${1}.job
  
done

echo "INFO: Resubmitted ${count} failed jobs for AmBeSource with ${1} as neutron absorber."

# Unset environment variables
unset RAT_SOURCE
unset DIR

