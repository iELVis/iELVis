#!/bin/sh

# iloc2fsurf.sh FsurfId iLocFolder
#
#The first argument is the patient's FreeSurfer Id (e.g., PT001)
#The second argument is the full path to the Patient's iLoc folder (e.g., /Users/Penfield/Desktop/PT001)
#
#This script uses FSL's flirt command to rigidly (i.e., 6 degrees of freedom mapping) transform the iLoc preimplant MRI so that it aligns with the FreeSurfer MRI. It then applies that same transformation to the postimplant scan and electrode coordinates. The result is an elec_recon subfolder in the patient's FreeSurfer folder that is then ready for iELVis's MATLAB plotting functions.
#Images of the two volumes overlayed are automatically generated so that you can inspect the quality of the coregistration.
#In the process the elec_recon subfolder in the patient's FreeSurfer folder is created along with the following files:
#    T1.nii.gz: The full head MRI
#    brainmask.nii.gz: The skull stripped MRI
#    postInPre.nii.gz: The post-implant CT coregistered to the pre-implant MRI
#
# Created by David Groppe on 2/11/15.
# Questions? Email: david.m.groppe@gmail.com
# Copyright 2019 Waveminer Analytics. All rights reserved.

usage='\nUSAGE:\n  iloc2fsurf.sh FsurfId iLocFolder\n\nEXAMPLE:\n iloc2fsurf.sh PT001 /Users/Penfield/Desktop/PT001\n'

if [ "$#" -ne 2 ]; then
 echo "Illegal number of parameters"
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

if [ ! -d  $2 ]; then
 echo
 echo "...Folder ${2} not found. Exit."
 echo
 exit 2
fi

elecReconPath=$SUBJECTS_DIR/$sub/elec_recon
echo 'Creating directory ' $elecReconPath
mkdir $elecReconPath

echo 'Copying iLoc electrode info file to elec_recon folder.'
cp $2/ElecInfo.tsv $elecReconPath/iLocElecInfo.tsv

echo 'Copying iLoc electrode pair file to elec_recon folder.'
cp $2/ElecPairs.tsv $elecReconPath/iLocElecPairs.tsv

echo 'Creating T1.nii.gz in elec_recon folder for coregistration.'
mriPath=$SUBJECTS_DIR/$sub/mri
mri_convert $mriPath/T1.mgz $elecReconPath/T1.nii.gz

echo 'Creating brainmask.nii.gz in elec_recon folder for potential visualization use later.'
mri_convert $mriPath/brainmask.mgz $elecReconPath/brainmask.nii.gz

echo 'Copying iLoc preimplant nii.gz file to elec_recon folder.'
cp $2/preimplant.nii.gz $elecReconPath/preimp_iLoc.nii.gz

echo 'Copying iLoc postimplant nii.gz file to elec_recon folder.'
cp $2/postinpre.nii.gz $elecReconPath/postimp_iLoc.nii.gz

echo 'Copying iLoc electrode info file to elec_recon folder.'
cp $2/ElecInfo.tsv $elecReconPath/iLocElecInfo.tsv

echo 'Registering ' $elecReconPath/preimp_iLoc.nii.gz ' to T1.nii.gz with a rigid (6 degrees of freedom) transformation that maximizes normalized correlation between the volumes. This takes awhile....'
flirt -in $elecReconPath/preimp_iLoc.nii.gz  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/iLocInFsurfPre.nii.gz -omat $elecReconPath/iLoc2fsurf.mat -interp trilinear -cost normcorr -dof 6 -searchcost normcorr -searchrx -180 180 -searchry -180 180 -searchrz -180 180

echo 'Now applying that same transformation to postimplant iLoc scan'
flirt -in $elecReconPath/postimp_iLoc.nii.gz  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/iLocInFsurfPost.nii.gz -interp trilinear -init $elecReconPath/iLoc2fsurf.mat -applyxfm
#flirt -in postopCT.nii.gz -ref T1.nii.gz -out ctINt1.nii.gz -interp trilinear -init ct2mriBBreg.mat -applyxfm

# Make directory to store coregistration images
mkdir -p $elecReconPath/PICS/COREG/

# Make images of MRI2MRI coregistration
slices $elecReconPath/iLocInFsurfPre.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/iLocInFsurfPre.nii.gz
# Make gifs of those images
slices $elecReconPath/iLocInFsurfPre.nii.gz  $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/iLocINfsurfT1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/iLocInFsurfPre.nii.gz  -o $elecReconPath/PICS/COREG/iLocINfsurfT1_2.gif

# Make images of CT2MRI coregistration
slices $elecReconPath/iLocInFsurfPost.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/iLocInFsurfPost.nii.gz
# Make gifs of those images
slices $elecReconPath/iLocInFsurfPost.nii.gz  $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/iLocINfsurfT1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/iLocInFsurfPost.nii.gz  -o $elecReconPath/PICS/COREG/iLocINfsurfT1_2.gif

echo 'Run the command below to interactively inspect the coregistration:'
echo "fslview ${elecReconPath}/T1.nii.gz ${elecReconPath}/iLocInFsurfPre.nii.gz"
