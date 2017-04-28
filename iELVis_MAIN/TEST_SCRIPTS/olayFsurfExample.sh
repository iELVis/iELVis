bbregister --s PT001 --mov /Users/davidgroppe/GIT/EpiSurg/iELVis/EXAMPLE_NII_FILES/fmriVolume.nii.gz --reg fmri2mri.dat --fslmat fmri2mri.mat --init-fsl --bold

tkregister2 --mov /Users/davidgroppe/GIT/EpiSurg/iELVis/EXAMPLE_NII_FILES/fmriVolume.nii.gz --reg fmri2mri.dat --surf

mri_vol2surf --src /Users/davidgroppe/GIT/EpiSurg/iELVis/EXAMPLE_NII_FILES/handMotorZVol.nii.gz --reg fmri2mri.dat --out handMotorLH.mgh --out_type mgh --hemi lh 

mri_vol2surf --src /Users/davidgroppe/GIT/EpiSurg/iELVis/EXAMPLE_NII_FILES/handMotorZVol.nii.gz --reg fmri2mri.dat --out handMotorRH.mgh --out_type mgh --hemi rh 


