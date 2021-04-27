#!/bin/sh

# gatherFsurfFiles.sh
# Collect the minimal FreeSurfer and BioimageSuite files for sharing. A minimal copy of the subject's
# FreeSurfer folder is created on the desktop and then zipped.
#
# Created by David Groppe on 4/2/15.
# Copyright 2015 __MyCompanyName__. All rights reserved.

usage='\nUSAGE:\n  gatherFsurfFiles.sh freesurferSubjectName \n\nEXAMPLE:\n gatherFsurfFiles.sh PT001\n'

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

# Create new fs folder on desktop
#newFolder=$1;
newFolder=~/Desktop/$1;
echo "Making directory ${newFolder} on the Desktop"
mkdir $newFolder
  
# Copy elec_recon stuff
mkdir $newFolder/elec_recon
cp $fsDir/elec_recon/*.nii.gz $newFolder/elec_recon/.
cp $fsDir/elec_recon/$1.mgrid $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.DURAL* $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.LEPTO* $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.PIAL* $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.INF* $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.CT* $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.POSTIMPLANT* $newFolder/elec_recon/.
cp $fsDir/elec_recon/*PostimpLoc*.txt $newFolder/elec_recon/.
cp $fsDir/elec_recon/*.electrodeNames $newFolder/elec_recon/.
cp $fsDir/elec_recon/*elec_recon.pdf $newFolder/elec_recon/.
cp $fsDir/elec_recon/00README.txt $newFolder/elec_recon/.

# Copy mri stuff
mkdir $newFolder/mri
mkdir $newFolder/mri/transforms
mkdir $newFolder/mri/orig
cp $fsDir/mri/brainmask.mgz $newFolder/mri/.
cp $fsDir/mri/aparc+aseg.mgz $newFolder/mri/.
cp $fsDir/mri/orig.mgz $newFolder/mri/.
cp $fsDir/mri/transforms/talairach.xfm $newFolder/mri/transforms/.
cp $fsDir/mri/orig/001.mgz $newFolder/mri/orig/.

# Copy build stamp
mkdir $newFolder/scripts
cp $fsDir/scripts/build-stamp.txt $newFolder/scripts/.

# Copy label stuff
mkdir $newFolder/label
cp $fsDir/label/*.aparc.a2009s.annot $newFolder/label/.
cp $fsDir/label/*.aparc.annot $newFolder/label/.
cp $fsDir/label/*Yeo2011*.mat $newFolder/label/.

# Copy surf stuff
mkdir $newFolder/surf
cp $fsDir/surf/*.curv $newFolder/surf/.
cp $fsDir/surf/*.inflated $newFolder/surf/.
cp $fsDir/surf/*.pial $newFolder/surf/.
cp $fsDir/surf/*.white $newFolder/surf/.
cp $fsDir/surf/*.sphere $newFolder/surf/.
cp $fsDir/surf/*.sphere.reg $newFolder/surf/.
cp $fsDir/surf/*.pial-outer-smoothed $newFolder/surf/.

# zip it!
zipName="${sub}fsurf.zip"
zip -r ~/Desktop/$zipName $newFolder



