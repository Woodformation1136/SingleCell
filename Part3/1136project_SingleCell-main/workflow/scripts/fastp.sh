#!/bin/bash

# Get arguments from command line ==============================================
while getopts "t:l:a:b:c:x:y:z:" args
do
    case "${args}" in
        t) THREADS=${OPTARG};;
        l) LOG_PREFIX=${OPTARG};;
        a) IN_READ1=${OPTARG};;
        b) IN_READ2=${OPTARG};;
        c) IN_READ=${OPTARG};;
        x) OUT_READ1=${OPTARG};;
        y) OUT_READ2=${OPTARG};;
        z) OUT_READ=${OPTARG};;
    esac
done


# Execute fastp ================================================================
set -euxo pipefail

usage() {
    echo "Usage: bash $0 [-t threads] <-s suffix> <-l log_prefix> \
    <-a input_read1 -b input_read2>|<-c input_read> \
    <-x output_read1 -y output_read2>|<-z output_read>" 1>&2 
}

exit_abnormal() {
    usage
    exit 1
}

fastp_paired_end () {
	fastp \
        --in1 ${1} \
		--in2 ${2} \
		--out1 ${3} \
		--out2 ${4} \
		--thread $THREADS \
		--detect_adapter_for_pe \
		--correction \
		--cut_front \
		--cut_tail \
		--disable_trim_poly_g
}

fastp_single_end () {
	fastp \
        --in1 ${1} \
		--out1 ${2} \
		--thread $THREADS \
		--cut_front \
		--cut_tail \
		--disable_trim_poly_g
}


#===============================================================================
# mkdir --parents ${LOG_DIR}/QC_documents

if [[ -n ${IN_READ1:-} ]] && [[ -n ${IN_READ2:-} ]]
then

    echo "===== Run fastp on "$IN_READ1" and "$IN_READ2" ====="
    date +%Y-%m-%d_%H:%M:%S

    echo "===== It's a PE reads. ====="
    fastp_paired_end $IN_READ1 $IN_READ2 $OUT_READ1 $OUT_READ2

    echo ""

    mv fastp.html ${LOG_PREFIX}.html
    mv fastp.json ${LOG_PREFIX}.json

elif [[ -n ${IN_READ:-} ]]
then

    echo "===== Run fastp on "$IN_READ" ====="
    date +%Y-%m-%d_%H:%M:%S

    echo "===== It's a SE reads. ====="
    fastp_single_end $IN_READ $OUT_READ

    echo ""

    mv fastp.html ${LOG_PREFIX}.html
    mv fastp.json ${LOG_PREFIX}.json

else

    exit_abnormal

fi

echo "====== Done ======"
date +%Y-%m-%d_%H:%M:%S
