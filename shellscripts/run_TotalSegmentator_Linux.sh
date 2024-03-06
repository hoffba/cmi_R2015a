#!/bin/bash
# Inputs:
# 	1: Image filename
# 	2: Saved segmentation filename
if [ $# == 2 ]; then
	echo ... Starting TS shell script

	module load python3.10-anaconda/2023.03
	conda init bash

	# First make sure python is set up to run TotalSegmentator
	CONDALIST=$(conda env list)
	if [[ $CONDALIST == *"TSenv"* ]]; then
		conda activate TSenv
	else
		echo ... Creating TSenv
		conda create -n TSenv
		conda activate TSenv
		conda install pytorch torchvision
		pip install totalsegmentator
		pip install cupy-cuda11x cucim
	fi
	USERNAME=$(whoami)
	TSFIND="$(find /home/"$USERNAME"/.conda/envs/TSenv -name "TotalSegmentator.py")"
	TSPATH="$(echo $TSFIND | tail -n1)"
	
	# Generate call to TotalSegmentator
	echo Input image: $1
	echo Save as: $2
	python $TSPATH -i $1 -o $2 --ml
	
elif [ $# -lt 2 ]; then
	echo "$0: Missing arguments"
	exit 1
else
	echo "$0: Too many arguments: $@"
	exit 1
fi