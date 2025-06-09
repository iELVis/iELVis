#!/bin/sh

# persyst2fsurf.sh FsurfId persystMinimalExportFolder
#
#The first argument is the patient's FreeSurfer Id (e.g., PT001)
#The second argument is the full path to the Patient's folder of Minimal-format Persyst files (e.g., /Users/Penfield/Desktop/PT001_Minimal). It should contain pre- and post-implant scans in approximate MNI space as well as ElecInfo.tsv, a tab separated text file of electrode coordinates, and ElecPairs.tsv, a tab separated text file that indicates which electrodes are neighbors and their color.
#
#This script uses FSL's flirt command to rigidly (i.e., 6 degrees of freedom mapping) transform Persyst's  preimplant MRI so that it aligns with the FreeSurfer MRI. It then applies that same transformation to the postimplant scan and electrode coordinates. The result is an elec_recon subfolder in the patient's FreeSurfer folder that is then ready for iELVis's MATLAB plotting functions.
#Images of the two volumes overlayed are automatically generated so that you can inspect the quality of the coregistration.
#In the process the elec_recon subfolder in the patient's FreeSurfer folder is created along with the following files (among others):
#    T1.nii.gz: The full head MRI
#    brainmask.nii.gz: The skull stripped MRI
#    postInPre.nii.gz: The post-implant CT coregistered to the pre-implant MRI
#
# Note also that FSL's img2imgcoord command is used to do the coordinate transformation

usage='\nUSAGE:\n  persyst2fsurf.sh FsurfId persystFolder\n\nEXAMPLE:\n persyst2fsurf.sh PT001 /Users/Penfield/Desktop/PT001_Minimal\n'

echo "persyst2fsurf.sh Version 2025-06-09 using mm coordinates"
echo " "

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

echo 'Copying Persyst electrode info file to elec_recon folder.'
cp $2/ElecInfo.tsv $elecReconPath/persystElecInfo.tsv

echo 'Copying Persyst electrode pair file to elec_recon folder.'
cp $2/ElecPairs.tsv $elecReconPath/persystElecPairs.tsv

echo 'Creating T1.nii.gz in elec_recon folder for coregistration.'
mriPath=$SUBJECTS_DIR/$sub/mri
mri_convert $mriPath/T1.mgz $elecReconPath/T1.nii.gz

echo 'Creating brainmask.nii.gz in elec_recon folder for potential visualization use later.'
mri_convert $mriPath/brainmask.mgz $elecReconPath/brainmask.nii.gz

echo 'Copying Persyst preimplant nii.gz file to elec_recon folder.'
cp $2/preimplant*.nii.gz $elecReconPath/preimp_persyst.nii.gz

echo 'Copying Persyst postimplant nii.gz file to elec_recon folder.'
cp $2/postInPre*.nii.gz $elecReconPath/postimp_persyst.nii.gz

echo 'Registering ' $elecReconPath/preimp_persyst.nii.gz ' to T1.nii.gz with a rigid (6 degrees of freedom) transformation that maximizes normalized correlation between the volumes. This takes awhile....'
flirt -in $elecReconPath/preimp_persyst.nii.gz  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/persystInFsurfPre.nii.gz -omat $elecReconPath/persyst2fsurf.mat -interp trilinear -cost normcorr -dof 6 -searchcost normcorr -searchrx -180 180 -searchry -180 180 -searchrz -180 180

echo 'Now applying that same transformation to postimplant Persyst scan'
flirt -in $elecReconPath/postimp_persyst.nii.gz  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/persystInFsurfPost.nii.gz -interp trilinear -init $elecReconPath/persyst2fsurf.mat -applyxfm
#flirt -in postopCT.nii.gz -ref T1.nii.gz -out ctINt1.nii.gz -interp trilinear -init ct2mriBBreg.mat -applyxfm

# Make directory to store coregistration images
mkdir -p $elecReconPath/PICS/COREG/

# Make images of MRI2MRI coregistration
slices $elecReconPath/persystInFsurfPre.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/persystInFsurfPre.nii.gz
# Make gifs of those images
slices $elecReconPath/persystInFsurfPre.nii.gz  $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/persystINfsurfT1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/persystInFsurfPre.nii.gz  -o $elecReconPath/PICS/COREG/persystINfsurfT1_2.gif

# Make images of CT2MRI coregistration
slices $elecReconPath/persystInFsurfPost.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/persystInFsurfPost.nii.gz
# Make gifs of those images
slices $elecReconPath/persystInFsurfPost.nii.gz  $elecReconPath/T1.nii.gz -o $elecReconPath/PICS/COREG/persystINfsurfT1_1.gif
slices $elecReconPath/T1.nii.gz  $elecReconPath/persystInFsurfPost.nii.gz  -o $elecReconPath/PICS/COREG/persystINfsurfT1_2.gif

echo 'Run the command below to interactively inspect the coregistration:'
echo "fsleyes ${elecReconPath}/T1.nii.gz ${elecReconPath}/persystInFsurfPre.nii.gz"

# Extract just the electrode coordinates
# Remove header and Extract columns 3-5
tail -n +2 $elecReconPath/persystElecInfo.tsv | cut -d'	' -f 3,4,5 > $elecReconPath/persystXyz.tsv

# Convert Persyst electrode coordinates. Note, img2imgcoord is an FSL function
img2imgcoord -src $elecReconPath/postimp_persyst.nii.gz -dest $elecReconPath/T1.nii.gz -xfm $elecReconPath/persyst2fsurf.mat $elecReconPath/persystXyz.tsv -mm > $elecReconPath/fsurfXyz.txt
