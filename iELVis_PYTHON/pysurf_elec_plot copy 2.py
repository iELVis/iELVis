"""
This script is based on the followying Pysurfer Source:
https://pysurfer.github.io/auto_examples/plot_basics.html#sphx-glr-auto-examples-plot-basics-py

You need to install Pysurfer (which uses Mayavi & Matplotlib for the visualization):

https://pysurfer.github.io/

NOTE: The plot is not 100% accurate. I think brain.add_foci assigns each electrode to the closest vertex on the
brain surface, which can distort the locations of some electrodes (e.g., the ones over sulci)

"""
import time
import pandas as pd
import numpy as np
from surfer import Brain
import os

# INPUT PARAMETERS (THIS IS WHAT YOU CHANGE)
sub='PT001'
BidsIeegRootDir='/Users/davidgroppe/Desktop/HandMotor/'
hemi = 'lh'

# Pysurfer/Freesurfer parameters
subject_id='sub-'+sub
subjects_dir=os.path.join(BidsIeegRootDir,'derivatives','iELVis')
surf = 'pial'


# Import electrode info
inFname=os.path.join(BidsIeegRootDir,subject_id,'ieeg',subject_id+'_ses-01_space-pial_electrodes.tsv')
xyzDf=pd.read_csv(inFname,sep='\t')
nElec=xyzDf.shape[0]
xyz=np.zeros((nElec,3))
# You need to subtract 1 from each coordinate to account for differences in Python and Matlab coordinates
# (since electrode coordinates were derived in Matlab)
for c in range(nElec):
    xyz[c,0]=xyzDf['x'][c]-1
    xyz[c,1]=xyzDf['y'][c]-1
    xyz[c,2]=xyzDf['z'][c]-1



"""
Call the Brain object constructor with these
parameters to initialize the visualization session.
"""
brain = Brain(subject_id, hemi, surf,subjects_dir=subjects_dir)


# Add electrodes
brain.add_foci(xyz, map_surface="white", color="gold",scale_factor=0.5)

input("Press Enter to continue...")

