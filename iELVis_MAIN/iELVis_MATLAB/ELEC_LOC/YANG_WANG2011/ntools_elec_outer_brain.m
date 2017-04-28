function ntools_elec_outer_brain(SubjectPath)
%function ntools_elec_outer_brain(SubjectPath)
%
% SubjectPath = Freesurfer subject folder

%  check the outer-brain lgi surface, if they are not there, create them

SubjectPath = [SubjectPath '/surf/'];
lgi_lh = exist([SubjectPath 'lh.pial-outer-smoothed'],'file');
lgi_rh = exist([SubjectPath 'rh.pial-outer-smoothed'],'file');

if lgi_lh==0
    % create the lh
    if exist([SubjectPath 'lh.pial'],'file');
        fprintf('Creating the left hemisphere outer smoothed brain surface \n');
        outer_smoothed('lh',SubjectPath);
    else
        fprintf('No lh.pial file detected, please finish the recon-all first \n');
    end
else
    fprintf('lh.pial-outer-smoothed file detected\n');
end

if lgi_rh==0
    % create the rh
    if exist([SubjectPath 'rh.pial'],'file');
        fprintf('Creating the right hemisphere outer smoothed brain surface \n');
        outer_smoothed('rh',SubjectPath);
    else
        fprintf('No rh.pial file detected, please finish the recon-all first \n');
    end
else
    fprintf('rh.pial-outer-smoothed file detected\n');
end

function outer_smoothed(sph,Sub_Path)

% this function is to create the pial-outer-smoothed surface by
% implementing first 4 steps of mris_compute_lgi

spath = [Sub_Path sph];

mris_fill = sprintf('mris_fill -c -r 1 %s.pial %s.pial.filled.mgz',spath,spath);
system(mris_fill);

% corr_oc = sprintf('use500; mri_convert %s.pial.filled.mgz -oc 0 0 0 %s.pial.filled.oc.mgz',spath,spath);
% system(corr_oc);
% make_outer_surface ([spath '.pial.filled.oc.mgz'],15,[spath '.pial-outer']);

make_outer_surface ([spath '.pial.filled.mgz'],15,[spath '.pial-outer']);

mris_extract = sprintf('mris_extract_main_component %s.pial-outer %s.pial-outer-main',spath,spath);
system(mris_extract);

mris_smooth = sprintf('mris_smooth -nw -n 30 %s.pial-outer-main %s.pial-outer-smoothed', spath, spath);
system(mris_smooth);


