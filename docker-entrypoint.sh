#!/bin/bash

function help {
cat <<EOF
Syntax: [options] [ANALYSIS]

networks    :: Predict tissue-specific functional networks
predictions :: Predict disease genes

Options:

-h | --help         :: this help

EOF
}

function predictions() {
    echo "Running disease prediction"
    if [[ "$#" -ne 2 ]]
    then
        echo "need disease, tissue network arguments"
        exit 1
    fi
    base_dir="/dq/data"
    dpred_dir="$base_dir/disease_prediction"

    disease_gold="$dpred_dir/gold_standard/${1}.label"
    tissue_net="$base_dir/network_integration/networks/${2}.dab"
    if [ ! -f $disease_gold ]
    then
        echo "disease not found"
        exit 1
    fi

    if [ ! -f $tissue_net ]
    then
        echo "tissue network not found"
        exit 1
    fi
    python /dq/src/predict_disease_genes.py -b /dq/bin/SVMperfer \
        -t $tissue_net \
        -l $disease_gold \
        -o /dq/outputs/${1}tmp \
        -p $dpred_dir/params/${1}.param
    mv /dq/outputs/${1}tmp* /dq/outputs/${1}.predictions
}

function networks() {
    echo "Running network integration"
    if [[ "$#" -ne 1 ]]
    then
        echo "need tissue name"
        exit 1
    fi
    echo $1
    base_dir="/dq/data/network_integration/"

    python /dq/src/network_integration.py -L -N -P \
      -B /dq/bin \
      -a $base_dir/gold_standard/global.dab \
      -A $1 \
      -z $base_dir/zeros.txt \
      -e $base_dir/genes.txt \
      -r $base_dir/alphas.txt \
      -d $base_dir/data_compendium/ \
      -s $base_dir/gold_standard/ \
      -W $base_dir/weights/ \
      -F -n \
      -w /dq/outputs/
    python /dq/src/network_integration.py --reprior 0.5 \
      -B /dq/bin \
      -a $base_dir/gold_standard/global.dab \
      -A $1 \
      -d $base_dir/data_compendium/ \
      -s $base_dir/gold_standard/ \
      -W $base_dir/weights/ \
      -w /dq/outputs/

    chmod 777 -R outputs/all

}

while [ "$1" ]
do case "$1" in
    -h | --help | help ) help
        exit
	    ;;
    * ) break
	;;
   esac
done


if [ -z $1 ] || [ $1 = 'networks' ]
then
  set -- "${@:2}"
  networks $@
fi
if [ -z $1 ] || [ $1 = 'predictions' ]
then
  set -- "${@:2}"
  predictions $@
fi
