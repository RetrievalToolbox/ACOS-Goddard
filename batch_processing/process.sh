#!/bin/bash

# Input environment variables (REQUIRED!)
# ORBIT: orbit string (don't forget the trailing character, if present!)
# SRCDIR: path to ACOS-Goddard repo
# SPECDIR: path to directory with spectroscopy files
# OUTDIR: path to directory in which results get stored


# We check if these variables exist, otherwise quit.



# Each of NJOB processes will be given a set of NCHUNK IDs at a time,
# so ./run_multiple will ingest NCHUNK IDs into ACOS-Goddard, which boots up
# once to process those 10 IDs. Then, the next chunk is processed, requiring
# another boot-up / call to ACOS-Goddard to process another NCHUNK IDs.


NJOBS=4 # How many jobs / processes to run in parallel
NLIST=8 # How many
NCHUNK=2 # How many chunks are processed by each parallel process in one go


ORBIT=${1}

# Find the L1b, L2Met and L2CPr files based on the orbit number
f_l1b=`ls ${SRCDIR}/data/*L1bSc*${ORBIT}*.h5`
N_l1b=`echo ${f_l1b} | wc -l`
if [ "$N_l1b" -eq 1 ]; then
    echo "Found L1b file at: ${f_l1b}"
else
    echo "None or more than one file found: ${f_l1b}"
    echo "Exiting"
    exit
fi

f_l2met=`ls ${SRCDIR}/data/*L2Met*${ORBIT}*.h5`
N_l2met=`echo ${f_l2met} | wc -l`
if [ "$N_l2met" -eq 1 ]; then
    echo "Found L2Met file at: ${f_l2met}"
else
    echo "None or more than one file found: ${f_l2met}"
    echo "Exiting"
    exit
fi

f_l2cpr=`ls ${SRCDIR}/data/*L2CPr*${ORBIT}*.h5`
N_l2cpr=`echo ${f_l2met} | wc -l`
if [ "$N_l2cpr" -eq 1 ]; then
    echo "Found L2CPr file at: ${f_l2cpr}"
else
    echo "None or more than one file found: ${f_l2cpr}"
    echo "Exiting"
    exit
fi

export L1B_FILE=${f_l1b}
export L2MET_FILE=${f_l2met}
export L2CPR_FILE=${f_l2cpr}
export SPECDIR=${SPECDIR}
export ORBIT=${ORBIT}
export OUTDIR=${OUTDIR}
export XRTM_PATH=${XRTM_PATH}


h5ls -d ${L1B_FILE}/SoundingGeometry/sounding_id \
    | grep -oP "\d{16}" \
    | head -n ${NLIST} \
    | parallel --keep-order -j ${NJOBS} -n ${NCHUNK} --progress ./run_multiple.sh
