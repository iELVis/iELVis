#!/bin/sh

# mri2mri.sh
# 
#The second argument is the nii.gz file of the postimplant MRI scan that you need to create with something like Matlab or FSL's mri_convert function.
#This script uses FSL's flirt command to rigidly (i.e., 6 degrees of freedom mapping) transform the CT scan so that it lines up with the preimplant MRI by maximizing the mutual information between the volumes.
#Images of the two volumes overlayed are automatically generated so that you can inspect the quality of the coregistration. 
#In the process the elec_recon subfolder in the patient's FreeSurfer folder is created along with the following files:
#    T1.nii.gz: The full head MRI
#    brainmask.nii.gz: The skull stripped MRI
#    postInPre.nii.gz: The post-implant CT coregistered to the pre-implant MRI
#
# Created by David Groppe on 2/11/15.
# Questions? Email: david.m.groppe@gmail.com
# Copyright 2015 __MyCompanyName__. All rights reserved.

usage='\nUSAGE:\n  mri2mri.sh freesurferSubjectName postimpMriNiiFile\n\nEXAMPLE:\n mri2mri.sh TWH014 /Users/dgroppe/Desktop/TWH_14_DICOMS/postimpMri.nii.gz\n'

if [ "$#" = 0 ]; then
 echo $usage
 exit 2
fi

sub=$1
fsDir=$SUBJECTS_DIR/$sub
if [ ! -d $fsDir ]; then
  echo
  echo "... ${fsDir} is not a directory."
  echo "...you have the wrong FreeSurfer codename for the subject or you have not yet run recon-all on this subject's MRI " 
  echo
  exit 2
fi

if [ ! -f  $2 ]; then
 echo
 echo "...File ${2} not found. Exit."
 echo
 exit 2
fi

elecReconPath=$SUBJECTS_DIR/$sub/elec_recon
echo 'Creating directory ' $elecReconPath
mkdir $elecReconPath

echo 'Creating T1.nii.gz in elec_recon folder for coregistration.'
mriPath=$SUBJECTS_DIR/$sub/mri
mri_convert $mriPath/T1.mgz $elecReconPath/T1.nii.gz

echo 'Creating brainmask.nii.gz in elec_recon folder for use with BioImageSuite later.'
mri_convert $mriPath/brainmask.mgz $elecReconPath/brainmask.nii.gz

echo 'Copying postimplant MRI nii.gz file to elec_recon folder.'
cp $2 $elecReconPath/postimpRaw.nii.gz

echo 'Registering ' $2 ' to T1.nii.gz with a rigid (6 degrees of freedom) transformation that maximizes normalized correlation between the volumes. This takes awhile....'
flirt -in $elecReconPath/postimpRaw.nii.gz  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/postInPre.nii.gz -omat $elecReconPath/post2preT1.mat -interp trilinear -cost normcorr -dof 6 -searchcost normcorr -searchrx -180 180 -searchry -180 180 -searchrz -180 180
# Make directory to store coregistration images
mkdir -p $elecReconPath/PICS/COREG/

# Make images of CT/MRI coregistration
slices $elecReconPath/postInPre.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/postInPre.nii.gz

# Make gifs of those images
slices $elecReconPath/postInPre.nii.gz $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/postINpreT1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/postInPre.nii.gz -o $elecReconPath/PICS/COREG/postINpreT1_2.gif

echo flirt -in $elecReconPath/postimpRaw.nii.gz  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/postInPre.nii.gz -omat $elecReconPath/post2preT1.mat -interp trilinear -cost normcorr -dof 6 -searchcost normcorr -searchrx -180 180 -searchry -180 180 -searchrz -180 180
echo 'Run the command below to interactively inspect the coregistration:'
echo "fslview ${elecReconPath}/T1.nii.gz ${elecReconPath}/postInPre.nii.gz"
