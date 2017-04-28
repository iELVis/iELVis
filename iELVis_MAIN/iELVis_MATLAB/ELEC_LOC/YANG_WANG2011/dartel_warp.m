% List of open inputs
% Segment: Data - cfg_files
% Segment: Tissue probability maps - cfg_files
% Initial Import: Output Directory - cfg_files
% Run DARTEL (existing Templates): Template - cfg_files
% Run DARTEL (existing Templates): Template - cfg_files
% Run DARTEL (existing Templates): Template - cfg_files
% Run DARTEL (existing Templates): Template - cfg_files
% Run DARTEL (existing Templates): Template - cfg_files
% Run DARTEL (existing Templates): Template - cfg_files
% Deformations: Image to base Id on - cfg_files
% Deformations: Apply to - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/home/lab/matlab/ntools_elec/dartel/dartel_warp_job_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(11, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Data - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Tissue probability maps - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Initial Import: Output Directory - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run DARTEL (existing Templates): Template - cfg_files
    inputs{5, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run DARTEL (existing Templates): Template - cfg_files
    inputs{6, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run DARTEL (existing Templates): Template - cfg_files
    inputs{7, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run DARTEL (existing Templates): Template - cfg_files
    inputs{8, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run DARTEL (existing Templates): Template - cfg_files
    inputs{9, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run DARTEL (existing Templates): Template - cfg_files
    inputs{10, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Image to base Id on - cfg_files
    inputs{11, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Apply to - cfg_files
end
spm('defaults', 'PET');
spm_jobman('serial', jobs, '', inputs{:});
