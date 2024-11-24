#!/usr/bin/env bash

folders=("yamanishi2008" "chembl" "wu2017_global")
fingerprints=("ecfp4") # 0.2
#fingerprints=("fcfp4") # 0.3
#fingerprints=("maccs") # 0.6

declare -A datasets
datasets["yamanishi2008"]="e gpcr ic nr"
datasets["chembl"]="Chembl28CCandD"
datasets["wu2017_global"]="global"

function run_loo()
{
    data_folder=$1
    data_set=$2
    fingerprint=$3
    
    echo "En: $data_folder, $data_set, $fingerprint"
    julia SimSpread.LOO.jl \
                        --dt data/${data_folder}/DT/${data_set}_DT.txt \
                        --dd data/${data_folder}/DD/${data_set}_DD.${fingerprint}_tanimoto.txt \
                        -o loocv-results/${data_set}-loocv/${fingerprint}/ \
                        --weighted
}

mkdir loocv-results
for folder in "${folders[@]}"; do
    for dataset in ${datasets[$folder]}; do
        mkdir loocv-results/${dataset}-loocv/
        for fingerprint in "${fingerprints[@]}"; do
            mkdir loocv-results/${dataset}-loocv/${fingerprint}
            run_loo $folder $dataset $fingerprint
        done
    done
done