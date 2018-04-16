#!/usr/bin/python

# takes as input a network + a disease label file
# runs svmperfer + makes predictions

import argparse
import subprocess

parser = argparse.ArgumentParser(description="given a tissue network and disease gene labels, makes predictions")

parser.add_argument('-t', '--tissue_network', required=True,
                    type=str,
                    help="tissue network dab")
parser.add_argument('-l', '--labels', required=True,
                    type=str,
                    help="disease gene labels (tab-delimited with gene <tab> label; 1 = pos, -1 = neg)")
parser.add_argument('-o', '--out', required=True,
                    type=str,
                    help="output file prefix")
parser.add_argument('-p', '--params_file',
                    type=str, required=True,
                    help="SVM params (see SVMperfer)")
parser.add_argument('-b', '--bin_svmperfer',
                    type=str, default="../bin/SVMperfer",
                    help="SVM binary location")

args = parser.parse_args()
subprocess.call([args.bin_svmperfer, '-i', args.tissue_network, '-l', args.labels, '-p', args.params_file, '-o', args.out, '-a'])
