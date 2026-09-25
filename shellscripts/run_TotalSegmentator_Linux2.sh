#!/bin/bash
# Inputs:
# 	1: Image filename
# 	2: Saved segmentation filename
if [ $# == 2 ]; then

	HOSTNAME=$(hostname)
	if [[ $HOSTNAME == "galban-ap-ps1a"* ]]; then
		echo ... Starting TS shell script on galban-ap-ps1a
		ml python-anaconda3/2023-04
	elif [[ $HOSTNAME == *".arc-ts."* ]]; then
		echo ... Starting TS shell script on Great Lakes
		ml python3.10-anaconda/2023.03
	else
		echo ... Unrecognized $HOSTNAME
		exit 1
	fi
		
	# Restart the shell
	conda init bash
	exec "$SHELL"

	# First make sure python is set up to run TotalSegmentator
	ENVPATH="/nfs/turbo/umms-cgalban/PyEnvironments/TSenv"
	if test -d $ENVPATH; then
		echo ... TSenv found
		conda activate $ENVPATH
	else
		echo ... Creating TSenv
		conda create -y -p $ENVPATH
		conda activate $ENVPATH
		conda install -y pytorch torchvision
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