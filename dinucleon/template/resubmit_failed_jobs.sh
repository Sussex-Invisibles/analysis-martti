#!/bin/bash
module load sge

# Check input arguments
if [ "$#" -ne 2 ]; then
  echo "ERROR: Need exactly 2 arguments (Decay mode, number of jobs)"
  return 1
fi

# Check for RAT installation
export RAT_SOURCE=/its/home/mn372/Software/SNOP/env_dinucleon.sh
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
export DIR=/its/home/mn372/Software/SNOP/work/analysis/dinucleon/${1}
if [ ! -d "$DIR" ]; then
  mkdir $DIR
fi
if [ ! -d "$DIR" ]; then
  echo "ERROR: Target directory does not exist: ${DIR}"
  return 4
fi
echo "Target directory is ${DIR}"

# Loop over bash logs
count=0
#for bashlog in ${DIR}/logs/bash/*; do
for index in `seq 1 ${2}`; do

  bashlog=${DIR}/logs/bash/*.${index}
  if [ -f ${bashlog} ] ; then

    # Skip jobs that completed successfully
    segfault=`tail ${bashlog} -n 2 | grep Segmentation | wc -l`
    last=`tail ${bashlog} -n 1`
    if [[ ${segfault} == 1 ]]; then
      echo "Segfault in job #${index}, will reprocess!"
    elif [[ ${last} == "Finished job script"* ]]; then
      continue
    fi
  fi

  # Reset/increase counters
  (( count++ ))
  found=0
  
  # Find bash log of failed job
  #index=`echo ${bashlog} | sed 's/.*\.//'`
  if [ -f ${bashlog} ] ; then 
    (( found++ ))
  else
    echo "WARNING: Could not find BASH log: ${bashlog}"
  fi
  
  # Find RAT log of failed job
  #line=`grep "Hostname: node" ${bashlog}`
  #node=`echo $line | cut -d " " -f 2 | sed 's/[^0-9]*//g'`
  #pid=`echo $line | cut -d " " -f 4`
  #ratlog=${DIR}/logs/rat_logs/rat.node${node}.${pid}.log
  ratlog=`grep ${1}_decay_${index}.root -r ${DIR}/logs/rat_logs/*.log | head -n 1 | cut -d ":" -f 1`
  if [ -f ${ratlog} ] ; then
    (( found++ ))
  else
    echo "WARNING: Could not find RAT log: ${ratlog}"
  fi
  
  # Find ROOT output file of failed job
  rootfile=${DIR}/output/${1}_decay_${index}.root
  if [ -f ${rootfile} ] ; then
    (( found++ ))
  else
    echo "WARNING: Could not find ROOT file: ${rootfile}"
  fi
  
  # Delete files ONLY if all files were found
  if [ $found -ne 3 ]; then
    echo "WARNING: Missing output files ($found) for job #${index}. Will not delete anything."
  else
    echo "Deleting files for job #${index}."
    rm ${bashlog}
    rm ${ratlog}
    #rm ${rootfile}
  fi
  
  # Resubmit jobs
  echo "qsub -jc mps.medium -q mps.q -l h_rt=48:00:00 -t ${index} ${DIR}/Dinucleon_${1}.job"
  #qsub -jc mps.medium -q mps.q -l h_rt=48:00:00 -t ${index} ${DIR}/Dinucleon_${1}.job
  
done

echo "INFO: Resubmitted ${count} failed jobs for Dinucleon Decay (${1} mode)."

# Unset environment variables
unset RAT_SOURCE
unset DIR

