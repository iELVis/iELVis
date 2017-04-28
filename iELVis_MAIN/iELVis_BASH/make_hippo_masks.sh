usage='\nUSAGE:\n  make_hippo_masks.sh freesurfer_subject\n\nEXAMPLE:\n make_hippo_masks.sh CaLi\n'

if [[ "$#" = 0 ]]; then
 echo -e $usage
 exit 2
fi

sub=$1
echo Making hippocampus masks for $1

#Convert aseg.mgz into nii file
echo mri_convert $SUBJECTS_DIR/$sub/mri/aseg.mgz  $SUBJECTS_DIR/$sub/mri/aseg.nii.gz
mri_convert $SUBJECTS_DIR/$sub/mri/aseg.mgz  $SUBJECTS_DIR/$sub/mri/aseg.nii.gz

#make hippo mask
echo 3dcalc -a $SUBJECTS_DIR/$sub/mri/aseg.nii.gz  -expr 'equals(a,17)+equals(a,53)' -prefix $SUBJECTS_DIR/$sub/elec_recon/hippo_mask.nii.gz
3dcalc -a $SUBJECTS_DIR/$sub/mri/aseg.nii.gz  -expr 'equals(a,17)+equals(a,53)' -prefix $SUBJECTS_DIR/$sub/elec_recon/hippo_mask.nii.gz
