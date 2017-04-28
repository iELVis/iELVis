usage='\nUSAGE:\n  make_hippo_masks.sh freesurfer_subject\n\nEXAMPLE:\n make_hippo_masks.sh CaLi\n'

if [[ "$#" = 0 ]]; then
 echo -e $usage
 exit 2
fi

sub=$1
echo Making filled.nii.gz file for $1

#Convert aseg.mgz into nii file
echo mri_convert $SUBJECTS_DIR/$sub/mri/filled.mgz  $SUBJECTS_DIR/$sub/mri/filled.nii.gz
mri_convert $SUBJECTS_DIR/$sub/mri/filled.mgz  $SUBJECTS_DIR/$sub/mri/filled.nii.gz
