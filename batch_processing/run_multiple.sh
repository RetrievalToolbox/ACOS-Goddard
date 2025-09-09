#!/bin/bash

# Needs the following ENVs set:
# SRCDIR: path to ACOS-Goddard repo
# SPECDIR: path to directory with spectroscopy files (TODO)
# L1B_FILE:
# L2MET:
# L2CPR:
# OUTDIR: path to output directory in which results are written into


# Check if the required ENV vars exist:
if [[ -z "${SRCDIR}" ]]; then echo "SRCDIR not set!"; exit; fi
if [[ -z "${SPECDIR}" ]]; then echo "SPECDIR not set!"; exit; fi
if [[ -z "${L1B_FILE}" ]]; then echo "L1B_FILE not set!"; exit; fi
if [[ -z "${L2MET_FILE}" ]]; then echo "L2MET_FILE not set!"; exit; fi
if [[ -z "${L2CPR_FILE}" ]]; then echo "L2CPR_FILE not set!"; exit; fi
if [[ -z "${OUTDIR}" ]]; then echo "OUTDIR not set!"; exit; fi


# Create output dir if it does not exist
if [ -d ${OUTDIR} ]; then
    echo "Output directory exists."
else
    echo "Creating output directory at ${OUTDIR}"
    mkdir -p ${OUTDIR}
fi

declare -a snid_process

# Check if any exist, they are not to be processed
for snid in "$@"; do

    if [ -f ${OUTDIR}/${snid}.h5 ]; then
        echo "File ${OUTDIR}/${snid}.h5 exists."
        continue
    fi
    if [ -f ${OUTDIR}/${snid}.h5.init ]; then
        echo "File ${OUTDIR}/${snid}.h5.init exists."
        continue
    fi
    if [ -f ${OUTDIR}/${snid}.h5.error ]; then
        echo "File ${OUTDIR}/${snid}.h5.error exists."
        continue
    fi

    # If they don't exist, we touch an init file and add it to the
    # list of those to be processed.

    touch ${OUTDIR}/${snid}.h5.init
    snid_process+=("${snid}")

done

# If this batch has no valid soundings to process, quit.
if [ ${#snid_process[@]} -eq 0 ]; then
    echo "No soundings in this batch to process."
    exit
fi

echo "This batch processes following IDs: ${snid_process[@]}"

echo "L1B: ${L1B_FILE}"
echo "L2Met: ${L2MET_FILE}"
echo "L2CPr: ${L2CPR_FILE}"
echo "Spectroscopy directory: ${SPECDIR}"
echo "Output directory: ${OUTDIR}"
echo "XRTM path: ${XRTM_PATH}"


cd ${SRCDIR}

MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=3 \
    julia --project=./ ${SRCDIR}/run.jl \
                  --solar_model ${SRCDIR}/example_data/l2_solar_model.h5 \
                  --L1b ${L1B_FILE} \
                  --L2Met ${L2MET_FILE} \
                  --L2CPr ${L2CPR_FILE} \
                  --o2_spec ${SPECDIR}/o2_v52.hdf \
                  --co2_spec ${SPECDIR}/co2_v52.hdf \
                  --h2o_spec ${SPECDIR}/h2o_v52.hdf \
                  --sounding_id ${snid_process[@]} \
                  --output ${OUTDIR} \
                  --spec 1,2,3 \
                  --polarized true \
                  --Lambertian false \
                  --aerosols true \
                  --retrieve_aerosols true \
                  --Nhigh 16 \
                  --LSI true \
                  --o2_scale 1.0048 \
                  --co2_scale_weak 0.994 \
                  --co2_scale_strong 0.998 \
                  --dsigma_scale 2.0 \
                  --gamma 50.0 \
                  --max_iterations 20


# Now loop again through all the procssed IDs and see which exist (at least executed successfully)
for snid in "${snid_process[@]}"; do

    echo "Checking ${snid}"
    if [ -f ${OUTDIR}/${snid}.h5 ]; then
        echo "Found ${snid} - success!"
        # SUCCESS!
        # remove init file move on
        rm ${OUTDIR}/${snid}.h5.init
    else
        echo "${snid} failed!"
        # BAD!
        # Remove init file and create error file
        rm ${OUTDIR}/${snid}.h5.init
        touch ${OUTDIR}/${snid}.h5.error
    fi

done
