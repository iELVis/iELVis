#!/bin/sh

# ct2mriBbreg.sh
#
# Registers a CT to an MRI using FreeSurfer's bbregister
#The second argument is the nii.gz file of the CT scan that you need to create with something like Matlab or FSL's mri_convert function.
#Images of the two volumes overlayed are automatically generated so that you can inspect the quality of the coregistration.
#In the process the elec_recon subfolder in the patient's FreeSurfer folder is created along with the following files:
#    T1.nii.gz: The full head MRI
#    brainmask.nii.gz: The skull stripped MRI
#    postInPre.nii.gz: The post-implant CT coregistered to the pre-implant MRI
#
# Created by David Groppe on 2/11/15.
# Questions? Email: david.m.groppe@gmail.com
# Copyright 2015 __MyCompanyName__. All rights reserved.

usage='\nUSAGE:\n  ct2mriBbreg.sh freesurferSubjectName ctNiiFile\n\nEXAMPLE:\n ct2mriBbreg.sh TWH014 /Users/dgroppe/Desktop/TWH_14_DICOMS/ct.nii.gz\n'

if [ "$#" = 0 ]; then
 echo -e $usage
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

echo 'Copying CT nii.gz file to elec_recon folder.'
cp $2 $elecReconPath/postimpRaw.nii.gz

bbregister --s $sub --mov postimpRaw.nii.gz --reg $elecReconPath/ct2mri.dat --fslmat $elecReconPath/ct2mri.mat --init-fsl --bold
flirt -in $elecReconPath/postimpRaw.nii.gz -ref $elecReconPath/T1.nii.gz -out $elecReconPath/postInPre.nii.gz -interp trilinear -init $elecReconPath/ct2mri.mat -applyxfm
# Make directory to store coregistration images
mkdir -p $elecReconPath/PICS/COREG/

# Make images of CT/MRI coregistration
slices $elecReconPath/postInPre.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/postInPre.nii.gz

# Make gifs of those images
slices $elecReconPath/postInPre.nii.gz $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/ctINt1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/postInPre.nii.gz -o $elecReconPath/PICS/COREG/ctINt1_2.gif

echo 'Run this for interactive GUI'
echo 'fslview ' $elecReconPath/T1.nii.gz ' ' $elecReconPath/postInPre.nii.gz
