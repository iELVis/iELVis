#!/bin/sh

# ct2mriManual.sh
#
# Registers a CT to an MRI using bbregister and a manual initilization stored in ct2mri.dat
#
#The second argument is the nii.gz file of the CT scan that you need to create with something like Matlab or FSL's mri_convert function.
#Images of the two volumes overlayed are automatically generated so that you can inspect the quality of the coregistration.
#This script expects an elec_recon subfolder in the patient's FreeSurfer folder with the following files:
#    T1.nii.gz: The full head MRI
#    ct2mri.mat: A manually adjusted affine transformation that aligns the CT scan with the preimplant MRI
#
#
# Created by David Groppe on 2/11/15.
# Questions? Email: david.m.groppe@gmail.com
# Copyright 2015 __MyCompanyName__. All rights reserved.

usage='\nUSAGE:\n  ct2mriManual.sh freesurferSubjectName ctNiiFile\n\nEXAMPLE:\n ct2mriBbreg.sh TWH014 /Users/dgroppe/Desktop/TWH_14_DICOMS/ct.nii.gz\n'


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

echo 'Copying CT nii.gz file to elec_recon folder.'
cp $2 $elecReconPath/postimpRaw.nii.gz

echo 'Copying ct2mri.dat to ct2mriManual.dat.'
cp $elecReconPath/ct2mri.dat $elecReconPath/ct2mriManual.dat

bbregister --s $sub --mov postimpRaw.nii.gz --reg $elecReconPath/ct2mri.dat --fslmat $elecReconPath/ct2mri.mat --init-reg $elecReconPath/ct2mriManual.dat  --bold
flirt -in ct2mriManual/postimpRaw.nii.gz -ref $elecReconPath/T1.nii.gz -out $elecReconPath/postInPre.nii.gz -interp trilinear -init $elecReconPath/ct2mri.mat -applyxfm

# Make images of CT/MRI coregistration
slices $elecReconPath/postInPre.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/postInPre.nii.gz

# Make gifs of those images
slices $elecReconPath/postInPre.nii.gz $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/ctINt1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/postInPre.nii.gz -o $elecReconPath/PICS/COREG/ctINt1_2.gif

echo 'Run this for interactive GUI'
echo 'fslview ' $elecReconPath '/T1.nii.gz' $elecReconPath '/postInPre.nii.gz'
