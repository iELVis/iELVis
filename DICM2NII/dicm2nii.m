function varargout = dicm2nii(src, dataFolder, varargin)
% DICM2NII converts dicom files into NIfTI or img/hdr files. 
% 
% DICM2NII(dcmSource, niiFolder, outFormat, MoCoOption, subjName)
% 
% The input arguments are all optional:
%  1. source file or folder. It can be a zip or tgz file, a folder containing
%     dicom files, or other convertible files. It can also contain wildcards
%     like 'run1_*' for all files start with 'run1_'.
%  2. folder to save result files.
%  3. output file format:
%      0 or 'nii'           for single nii uncompressed.
%      1 or 'nii.gz'        for single nii compressed (default).
%      2 or 'hdr'           for hdr/img pair uncompressed.
%      3 or 'hdr.gz'        for hdr/img pair compressed.
%      4 or '3D.nii'        for 3D nii uncompressed (SPM12).
%      5 or '3D.nii.gz'     for 3D nii compressed.
%      6 or '3D.hdr'        for 3D hdr/img pair uncompressed (SPM8).
%      7 or '3D.hdr.gz'     for 3D hdr/img pair compressed.
%  4. MoCo series options:
%      0 create files for both original and MoCo series.
%      1 ignore MoCo series if both present (default).
%      2 ignore original series if both present.
%     Note that if only one of the two series is present, it will be created
%     always.
%  5. subject name for the data. The code can do only one subject each time.
%     If you have a folder or zip/tgz file for multiple subjects (not
%     recommended), you can specify which subject to convert with the 5th input.
%     The name should be the Patient Name entered at scanner console for in
%     format of LastName^FirstName, such as Smith^John, or simply the Last name
%     if no First name was provided to the console. If you already have mixed
%     subject data in a folder, you can let the code return the unconverted
%     subject names in the second output argument, and then provide each of
%     subjects as the 5th input to convert data into specific subject folder
%     (see example below).
%     A better method to deal with mixed subject data is to use SORT_DICM to
%     sort files into different subject folder, then convert one by one.
% 
% The optional output are converted (1st) and unconverted (2nd) PatientName.
% 
% Typical examples:
%  dicm2nii; % bring up user interface if there is no input argument
%  dicm2nii('D:/myProj/zip/subj1.zip', 'D:/myProj/subj1/data'); % zip file
%  dicm2nii('D:/myProj/subj1/dicom/', 'D:/myProj/subj1/data'); % folder
% 
% Less useful examples:
%  dicm2nii('D:/myProj/dicom/', 'D:/myProj/subj2/data', [], [], 'subj2');
%  dicm2nii('D:/myProj/dicom/run2*', 'D:/myProj/subj/data');
%  dicm2nii('D:/dicom/', 'D:/data', '3D.nii'); % SPM style files
% 
% Example to deal with multiple subjects:
%  [converted, otherSubj] = dicm2nii('D:/myProj/dicom/', 'D:/myProj/subj1/data');
%  movefile('D:/myProj/subj1', ['D:/myProj/' converted]); % rename folder
%  In case of multiple subjects, above command will convert one of the subjects,
%  and return the unconverted subjects in second output. You can convert other
%  subjects one by one using script like this:
%     for i = 1:length(otherSubj)
%         dataDir = ['D:/myProj/' otherSubj{i} '/data'];
%         dicm2nii('D:/myProj/dicom/', dataDir, [], [], otherSubj{i});
%     end
% 
% If there is no input, or any of the first two inputs is empty, a graphic user
% interface will appear.
% 
% If the first input is a zip/tgz file, such as those downloaded from dicom
% server, DICM2NII will extract files into a temp folder, create NIfTI files
% into the data folder, and then delete the temp folder. For this reason, it is
% better to keep the compressed file as backup.
% 
% If a folder is the data source, DICM2NII will convert all files in the folder
% and its subfolders (you don't need to sort files for different series).
% 
% Please note that, if a file in the middle of a series is missing, the series
% will be skipped without converting, and a warning message will be shown as red
% text in Command Windows, and saved into a text file the data folder.
% 
% Slice timing information, if available, is stored in nii header, such as
% slice_code and slice_duration. But the simple way may be to use the field
% SliceTiming in dcmHeaders.mat. That timing is actually those numbers for FSL
% when using custom slice timing. This is the universal method to specify any
% kind of slice order, and for now, is the only way which works for multiband.
% Slice order is one of the most confusing parameters, and it is recommended to
% use this method to avoid mistake. To convert this timing into slice order for
% SPM: [~, spm_order] = sort(s.SliceTiming, 'descend');
% 
% If there is DTI series, bval and bvec files will be generated for FSL etc. A
% Matlab data file, dcmHeaders.mat, is always saved into the data folder. This
% file contains dicom header from the first file for created series and some
% information from last file in field LastFile. For DTI series, B_value and
% DiffusionGradientDirection for all directions are saved into the mat file. For
% MoCo series, motion parameters, RBMoCoTrans and RBMoCoRot, are also saved.
% 
% Some information, such as TE, phase encoding direction and effective dwell
% time are stored in descrip of nii header. These are useful for fieldmap B0
% unwarp correction. Acquisition start time and date are also stored, and this
% may be useful if one wants to align the functional data to some physiological
% recording, like pulse, respiration or ECG.
% 
% The output file names adopt SeriesDescription or ProtocolName of each series
% used on scanner console. If both original and MoCo series are created, '_MoCo'
% will be appended for MoCo series. For phase image, such as those from field
% map, '_phase' will be appended to the name. In case of name conflict,
% SeriesNumber, such as '_s005', will be appended to make file names unique. It
% is suggested to use short and descriptive SeriesDescription on the scanner
% console, and use names containing only letters, numbers and underscores.
% 
% For SPM 3D files, the file names will have volume index in format of '_00001'
% appended to above name.
% 
% Some of the parameters may not be available for all vendors. For example,
% there is no slice order information in Philips data. Please report any bug to
% xiangrui.li@gmail.com or at
% http://www.mathworks.com/matlabcentral/fileexchange/42997
%
% DG modified this so that you can simply enter:
% dicm2nii
% and use GUIs to do the rest

% Thanks to:
% Jimmy Shen's Tools for NIfTI and ANALYZE image,
% Chris Rorden's dcm2nii pascal source code,
% Przemyslaw Baranski for direction cosine matrix to quaternions. 

% History (yymmdd):
% 130512 Publish to CCBBI users (Xiangrui Li).
% 130513 Convert img from uint16 to int16 if range allows;
%        Support output file format of img/hdr/mat.
% 130515 Change creation order to acquisition order (more natural).
%        If MoCo series is included, append _MoCo in file names.
% 130516 Use SpacingBetweenSlices, if exists, for SliceThickness. 
% 130518 Use NumberOfImagesInMosaic in CSA header (work for some old data).
% 130604 Add scl_inter/scl_slope and special naming for fieldmap.
% 130614 Work out the way to get EffectiveEchoSpacing for B0 unwarp.
% 130616 Add needed dicom field check, so it won't err later.
% 130618 Reorient if non-mosaic or slice_dim is still 3 and no slice flip.
% 130619 Simplify DERIVED series detection. No '_mag' in fieldmap name.
% 130629 Improve the method to get phase direction;
%        Permute img dim1&2 (no -90 rotation) & simplify xform accordingly.
% 130711 Make MoCoOption smarter: create nii if only 1 of 2 series exists.
% 130712 Remove 5th input (allHeader). Save memory by using partial header.
% 130712 Bug fix: dim_info with reorient. No problem since no EPI reorient.
% 130715 Use 2 slices for xform. No slice flip needed except revNum mosaic.
% 130716 Take care of lower/upper cases for output file names;
%        Apply scl_slope and inter to img if range allows and no rounding;
%        Save motion parameters, if any, into dcmHeader.mat.
% 130722 Ugly fix for isMos, so it works for '2004A 4VA25A' phase data;
%        Store dTE instead of TE if two TE are used, such as fieldmap.
% 130724 Add two more ways for dwell time, useful for '2004A 4VA25A' dicom.
% 130801 Can't use DERIVED since MoCoSeries may be labeled as DERIVED.
% 130807 Check PixelSpacing consistency for a series;
%        Prepare to publish to Matlab Central.
% 130809 Add 5th input for subjName, so one can choose a subject.
% 130813 Store ImageComments, if exists and is meaningful, into aux_file.
% 130818 Expand source to dicom file(s) and wildcards like run1*.dcm.
%        Update fields in dcmHeader.mat, rather than overwriting the file.
%        Include save_nii etc in the code for easy distribution.
% 130821 Bug fix for cellstr input as dicom source.
%        Change file name from dcm2nii.m to reduce confusion from MRICron.
%        GUI implemented into the single file.
% 130823 Remove dependency on Image Processing Toolbox.
% 130826 Bug fix for '*' src input. Minor improvement for dicm_hdr.
% 130827 Try and suggest to use pigz for compression (thanks Chris R.).
% 130905 Avoid the missing-field error for DTI data with 2 excitations.
%        Protect GUI from command line plotting.
% 130912 Use lDelayInTR for slice_dur, possibly useful for old data.
% 130916 Store B_matrix for DTI image, if exists.
% 130919 Make the code work for GE and Philips dicom at Chris R website.
% 130922 Remove dependence on normc from nnet toolbox (thank Zhiwei);
%        Prove no slice order info in Philips, at least for Intera 10.4.1.
% 130923 Make the code work for Philips PAR/REC pair files.
% 130926 Take care of non-mosaic DTI for Siemens (img/bval/bvec);
% 130930 Use verify_slice_dir subfun to get slice_dir even for a single file.
% 131001 dicm_hdr can deal with VR of SQ. This slows down it a little.
% 131002 Avoid fullfile for cellstr input (not supported in old ver matlab).
% 131006 Tweak dicm_hdr for multiframe dicom (some bug fixes);
%        First working version for multiframe (tested with Philips dicom).
% 131009 Put dicm_hdr, dicm_img, dicm_dict outside this file;
%        dicm_hdr can read implicit VR, and is faster with single fread;
%        Fix problem in gzipOS when folder name contains space.
% 131020 Make TR & ProtocolName non-mandatory; Set cal_min & cal_max.
% 131021 Check SamplesPerPixel, skip run if it is 1+.
% 131021 Implement conversion for AFNI HEAD/BRIK.
% 131024 Bug fix for dealing with current folder as src folder.
% 131029 Bug fix: Siemens, 2D, non-mosaic, rev-num slices were flipped.
% 131105 DTI parameters: field names more consistent; read DTI flds in
%        save_dti_para for GE/Philips (make others faster); convert Philips
%        bvec from deg into vector (need to be verified).
% 131114 Treak for multiframe dicm_hdr: MUCH faster by using only 1,2,n frames;
%        Bug fix for Philips multiframe DTI parameters;
%        Split multiframe Philips B0 map into mag and phase nii.
% 131117 Make the order of phase/mag image in Philips B0 map irrelevant.
% 131219 Write warning message to a file in data folder (Gui's suggestion).
% 140120 Bug fix in save_dti_para due to missing Manufacturer (Thank Paul).
% 140121 Allow missing instance at beginning of a series.
% 140123 save_nii: bug fix for gzip.m detection, take care of ~ as home dir.
% 140206 bug fix: MoCo detetion bug introduced by removing empty cell earlier.
% 140223 add missing-file check for Philips data by slice locations.
% 140312 use slice timing to set slice_code for both GE and Siemens.
%        Interleaved order was wrong for GE data with even number of slices. 
% 140317 Use MosaicRefAcqTimes from last vol for multiband (thank Chris).
%        Don't re-orient fieldmap, so make FSL happy in case of non_axial. 
%        Ugly fix for wrong dicom item VR 'OB': Avoid using main header 
%        in csa_header(), convert DTI parameters to correct type. There may
%        be other wrong parameters we don't realize. 
% 140319 Store SliceTiming field in dcmHeaders.mat for FSL custom slice timing.
%        Re-orient even if flipping slices for 2D MRAcquisitionType.
% 140324 Not set cal_min, cal_max anymore.
% 140327 Return unconverted subject names in 2nd output.
% 140401 Always flip image so phase dir is correct.
% 140409 Store nii extension (not enabled due to nifti ext issue).
% 140501 Fix for GE: use LocationsInAcquisition to replace ImagesInAcquisition;
%            isDTI=DiffusionDirection>0; Gradient already in image reference.
% 140505 Always re-orient DTI. bvec fix for GE DTI (thx Chris).
% 140506 Remove last DTI vol if it is computed ADC (as dcm2niix);
%        Use SeriesDescription to replace ProtocolName for file name;
%        Improved dim_info and phase direction.
% 140512 Decode GE ProtocolDataBlock for phase direction;
%        strtrim SeriesDescription for nii file name.
% 140513 change stored phase direction to image space for FSL unwarp;
%        Simplify code for dim_info.
% 140516 Switch back to ProtocolName for SIEMENS to take care of MOCO series;
%        Detect Philips Dim3IsVolume (for multi files) during dicom check; 
%        Work for GE interleaved slices even if InstanceNumber is in time order;
%        Do ImagePositionPatient check for all vendors;
%        Simplify code for save_dti_para.
% 140517 Store img with first dim flipped, to take care of DTI bvec problems. 
% 140522 Use SliceNormalVector for mosaic slice_dir, so no worry for revNumb;
%        Bug fix for interleaved descending slice_code.
% 140525 xform sliceCenter to SliceLocation in verify_slice_dir. 
% 140526 Take care of non-unique ixyz. 
% 140608 Bug fix for GE interleaved slices;
%        Take care all ixyz, put verify_slice_dir into xform_mat.
% 140610 Compute readout time for DTI, rather than dwell time.
% 140621 Support tgz file as data source.
% 140716 Bug fix due to empty src for GUI subject option.
% 140808 Simplify mosaic detection, and remove isMosaic.
% 140816 Simplify DTI detection.
% 140911 Minor fix for Siemens ProtocolName for error message.
% 141016 Remember GUI settings from last conversion;
%        Make multi-subject error message friendly.
% 141021 Show percent progress for validating dicom files.
% 141023 Get LocationsInAcquisition for GE multiframe dicom.
% 141024 Use unique ImagePositionPatient to determine LocationsInAcquisition.
% 141028 Use matlabpool if available and worthy.
% 141125 Store NumberOfTemporalPositions in dicom header.
% 141128 Minor tweaks for Octave 3.8.1 command line (GUI not working).
% 141216 Use ImagePositionPatient to derive SliceThickness if possible.
% 141217 Override LocationsInAcquisition with computed nSL (thx Luigi);
%        Check RescaleIntercept and RescaleSlope consistency.
% 141218 Allow 1e-4 diff for ImagePositionPatient of same slice location.
% 141223 multiFrameFields: return earlier if only single frame (thx Sander);
%        No re-orient for single slice (otherwise problem for mricron to read).
% 141224 mos2vol: use nSL loop (faster unless many slices).
% 141229 Save nii ext (ecode=40) if FSL is detected & it is not 5.0.5/5.0.6.
% 141230 nojvm: no matlabpool; no dicm_hdr progress due to '\b' issue for WIN.
% 150109 dicm_img(s, 0) to follow the update for dicm_img.
% 150112 Use nii_tool.m, remove make_nii, save_nii etc from this file.
% 150115 Allow SamplesPerPixel>1, but likely not very useful.
% 150117 Store seq name in intent_name.
% 150119 Add phase img detection for Philips (still need it for GE).
% 150120 No file skip by EchoTime: keep all data by using EchoNumber.
% 150209 Add more output format for SPM style: 3D output;
%        GUI includes SPM 3D, separates GZ option. 
% 150211 No missing file check for all vendors, relying on ImagePosition check;
%        csa_header() relies on dicm_hdr decoding (avoid error on old data);
%        Deal with dim3-RGB and dim4-frames due to dicm_img.m update.
% 150222 Remove useless, mis-used TriggerTime for partial hdr; also B_matrix.
% 150302 No hardcoded sign change for DTI bvec, except for GE;
%        set_nii_header: do flip only once after permute;
% 150303 Bug fix for phPos: result was right by lucky mistake;
%        Progress shows nii dim, more informative than number of files.
% End of history. Don't edit this line!

% DG added below
if nargin<1,
    [filename, src] = uigetfile('*.dcm', 'Pick a MATLAB code file');
    if isequal(filename,0) || isequal(src,0)
        disp('User pressed cancel')
    else
        disp(['User selected ', fullfile(src, filename)])
    end
    dataFolder='/Users/dgroppe/Desktop';
end

if nargin>1 && ischar(src) && strcmp(src, 'dicm2nii_gui_cb')
    dicm2nii_gui(dataFolder); % mis-use first two input for GUI
    varargout = {'' ''};
    return;
end

%% Deal with output format first, and error out if invalid
if nargin<3 || isempty(varargin{1}), fmt = 1; % default .nii.gz
else fmt = varargin{1};
end

if (isnumeric(fmt) && any(fmt==[0 1 4 5])) || ...
      (ischar(fmt) && ~isempty(regexpi(fmt, 'nii')))
    ext = '.nii';
elseif (isnumeric(fmt) && any(fmt==[2 3 6 7])) || (ischar(fmt) && ...
        (~isempty(regexpi(fmt, 'hdr')) || ~isempty(regexpi(fmt, 'img'))))
    ext = '.img';
else
    error(' Invalid output file format (the 3rd input).');
end

if (isnumeric(fmt) && mod(fmt,2)) || ...
        (ischar(fmt) && ~isempty(regexpi(fmt, '.gz')))
    ext = [ext '.gz']; % gzip file
end

rst3D = (isnumeric(fmt) && fmt>3) || ...
    (ischar(fmt) && ~isempty(regexpi(fmt, '3D')));

%% Deal with MoCo option
if nargin<4 || isempty(varargin{2})
    MoCo = 1; % by default, use original series if both present 
else
    MoCo = varargin{2};
    if ~any(MoCo==0:2)
        error(' Invalid MoCoOption. The 4th input must be 0, 1 or 2.');
    end
end

%% Deal with 5th input: we do one subject once
if nargin<5 || isempty(varargin{3})
    subjProvided = false; subj = '';
else 
    subjProvided = true; subj = varargin{3};
    if ~ischar(subj), error(' Invalid subject name.');end
end

%% Deal with data source
varargout = {};
unzip_cmd = '';
if isempty(src) || isempty(dataFolder)
% ORIG if nargin<1 || isempty(src) || (nargin<2 || isempty(dataFolder))
    create_gui; % show GUI if input is not enough
    return;
end

if isnumeric(src)
    error('Invalid dicom source.');    
elseif iscellstr(src) % multiple files
    dcmFolder = folderFromFile(src{1});
    n = length(src);
    fnames = src;
    for i = 1:n
        foo = dir(src{i});
        if isempty(foo), error('%s does not exist.', src{i}); end
        fnames{i} = fullfile(dcmFolder, foo.name); 
    end
elseif ~exist(src, 'file') % like input: run1*.dcm
    fnames = dir(src);
    if isempty(fnames), error('%s does not exist.', src); end
    fnames([fnames.isdir]) = [];
    dcmFolder = folderFromFile(src);
    fnames = strcat(dcmFolder, filesep, {fnames.name});    
elseif isdir(src) % folder
    dcmFolder = src;
elseif ischar(src) % 1 dicom or zip/tgz file
    dcmFolder = folderFromFile(src);
    unzip_cmd = compress_func(src);
    if isempty(unzip_cmd)
        fnames = dir(src);
        fnames = strcat(dcmFolder, filesep, {fnames.name});
    end
else 
    error('Unknown dicom source.');
end
dcmFolder = fullfile(getfield(what(dcmFolder), 'path'));

%% Deal with dataFolder
if nargin<2 || isempty(dataFolder)
    dataFolder = uigetdir(dcmFolder, 'Select a folder to save data files');
    if dataFolder==0, return; end
end
if ~isdir(dataFolder), mkdir(dataFolder); end
dataFolder = fullfile([getfield(what(dataFolder), 'path') filesep]);
global dcm2nii_errFileName;
if isempty(dcm2nii_errFileName) % show once for a session
    disp('Xiangrui Li''s dicm2nii (feedback to xiangrui.li@gmail.com)');
end
dcm2nii_errFileName = [dataFolder 'dicm2nii_warningMsg.txt'];

%% Unzip if compressed file is the source
tic;
if ~isempty(unzip_cmd)
    [~, fname, ext1] = fileparts(src);
    dcmFolder = sprintf('%stmpDcm%s/', dataFolder, fname);
    if ~isdir(dcmFolder), mkdir(dcmFolder); end
    disp(['Extracting files from ' fname ext1 ' ...']);

    if strcmp(unzip_cmd, 'unzip')
        cmd = sprintf('unzip -qq -o %s -d %s', src, dcmFolder);
        err = system(cmd); % first try system unzip
        if err, unzip(src, dcmFolder); end % Matlab's unzip is too slow
    elseif strcmp(unzip_cmd, 'untar')
        if isempty(which('untar')), error('No untar found in matlab path.'); end
        untar(src, dcmFolder);
    end
    drawnow;
end 

%% Get all file names including those in subfolders, if not specified
if ~exist('fnames', 'var')
    dirs = genpath(dcmFolder);
    dirs = textscan(dirs, '%s', 'Delimiter', pathsep);
    dirs = dirs{1}; % cell str
    fnames = {};
    for i = 1:length(dirs)
        curFolder = [dirs{i} filesep];
        foo = dir(curFolder); % all files and folders
        foo([foo.isdir]) = []; % remove folders
        foo = strcat(curFolder, {foo.name});
        fnames = [fnames foo]; %#ok<*AGROW>
    end
end
nFile = length(fnames);
if nFile<1, error(' No files found in the data source.'); end

%% Get Manufacturer
dict = dicm_dict('', 'Manufacturer');
vendor = '';
for i = unique([1 ceil(nFile*[0.2 0.5 0.8 1])]) % try up to 5 files
    s = dicm_hdr(fnames{i}, dict);
    if isempty(s), continue; end
    vendor = strtok(s.Manufacturer); % take 1st word only
    break;
end

%% Check each file, store header in cell array hh
% first 6 fields are must for 1st round, next 3 are must for later check
flds = {'InstanceNumber' 'SeriesNumber' 'ImageType' 'Columns' 'Rows' ...
	'BitsAllocated' 'PixelSpacing' 'ImageOrientationPatient' ...
    'ImagePositionPatient' 'PixelRepresentation' 'BitsStored' 'HighBit' ...
    'SamplesPerPixel' 'PlanarConfiguration' 'SeriesDescription' ...
    'EchoNumber' 'PatientName' 'PatientID' 'NumberOfFrames' ...
    'B_value' 'DiffusionGradientDirection' ...
    'RTIA_timer' 'RBMoCoTrans' 'RBMoCoRot' 'RescaleIntercept' 'RescaleSlope' };
dict = dicm_dict(vendor, flds); % get partial dict, suppose one vendor only
% Following for Philips only: B_value etc may be duplicated in later tags
ind = find(strcmp(dict.name, 'RescaleSlope'), 1, 'last') + 1;
dict.name(ind:end) = []; dict.tag(ind:end) = []; dict.vr(ind:end) = []; 

junk = {'\MEAN' '\DUMMY IMAGE' '\TTEST' '\FMRI\DESIGN' ... % GLM
        '\DIFFUSION\ADC\' '\DIFFUSION\FA\' '\DIFFUSION\TRACEW\'}; % DTI

% read header for all files, use matlabpool if available and worthy
hh = cell(1, nFile); errStr = cell(1, nFile);
doParal = usejava('jvm') && nFile>1000;
if doParal
    try 
        newOpen = matlabpool('size')<1;
        if newOpen
            try
                matlabpool;
            catch me
                fprintf(2, '%s\n', me.message);
                doParal = false;
            end
        end
    catch
        doParal = false;
    end
end

fprintf('Validating %g files (%s) ...\n', nFile, vendor);
if doParal
    try
        tObj = timerfindall('TimerFcn', @closeMatlabpool); 
        if newOpen && isempty(tObj)
            tObj = timer('StartDelay', 3600, ... % close matlabpool
                'TimerFcn', @closeMatlabpool, 'ObjectVisibility', 'off');
        elseif ~isempty(tObj) && strcmpi(tObj.Running, 'on')
            stop(tObj);
        end
        start(tObj);
    end
    parfor k = 1:nFile, [hh{k}, errStr{k}] = dicm_hdr(fnames{k}, dict); end
else
    for k = 1:nFile, disp(['k=' num2str(k)]); [hh{k}, errStr{k}] = dicm_hdr(fnames{k}, dict); end
end

%% sort headers into cell h by SeriesNumber, EchoNumber and InstanceNumber
h = {}; % in case of no dicom files at all
subj_skip = {};
errInfo = '';
for k = 1:nFile
    s = hh{k};
    if isempty(s) || any(~isfield(s, flds(1:6))) || isType(s, junk)
        if ~isempty(errStr{k})
            errInfo = sprintf('%s\n%s\n', errInfo, errStr{k});
        end
        continue; % skip the file
    end
    subj1 = tryGetField(s, 'PatientName');
    if isempty(subj1), subj1 = tryGetField(s, 'PatientID', 'unknown'); end
       
    % if not for single subject, do the first only
    if isempty(subj)
        subj = subj1; % store it for later check
    elseif ~strcmpi(subj, subj1)
        if ~any(strcmp(subj_skip, subj1)), subj_skip{end+1} = subj1; end
        continue;
    end

    i = tryGetField(s, 'EchoNumber', 1);  if i<1, i = 1; end
    j = s.InstanceNumber; if j<1, j = 1; end
    h{s.SeriesNumber}{i}{j} = s; % store partial header
end
clear hh errStr;

if nargout>0, varargout{1} = subj; end % return converted subject ID
if nargout>1, varargout{2} = unique(subj_skip); end % unconverted subject ID
if ~subjProvided && ~isempty(subj_skip)
    subj_skip = sprintf('%s, ', subj_skip{:});
    errorLog(['Files for following subjects are not converted: ' ...
        subj_skip(1:end-2)]);
end

%% Check headers: remove file-missing and dim-inconsistent series
nRun = length(h);
if nRun<1
    errorLog(sprintf('No valid files found for %s:\n%s.', subj, errInfo)); 
    return;
end
keep = true(1, nRun); % true for useful runs
isMoCo = false(1, nRun); % deal moco together later
fldsCk = {'ImageOrientationPatient' 'NumberOfFrames' 'Columns' 'Rows' ...
          'PixelSpacing' 'RescaleIntercept' 'RescaleSlope' 'SamplesPerPixel'};
for i = 1:nRun
    if isempty(h{i}), keep(i) = 0; continue; end % must be here due to MoCo
    h{i} = [h{i}{:}]; % concatenate different EchoNumber
    ind = cellfun(@isempty, h{i});
    h{i}(ind) = []; % remove all empty Instance for all vendors
    
    s = h{i}{1};
    series = sprintf('Subject %s, %s (Series %g)', subj, ProtocolName(s), s.SeriesNumber);
    if ~isfield(s, 'LastFile') % avoid re-read for PAR/AFNI file
        s = dicm_hdr(s.Filename); % full header for 1st file
    end
    s = multiFrameFields(s); % no-op if non multi-frame
    if ~isfield(s, 'Manufacturer'), s.Manufacturer = 'Unknown'; end
    h{i}{1} = s; % update record for full hdr or multiframe
    if any(~isfield(s, flds(7:9))), keep(i) = 0; continue; end
    isMoCo(i) = isType(s, '\MOCO\');
    
    % check consistency in 'fldsCk'
    nFile = length(h{i});
    if nFile<2, continue; end
    for k = 1:length(fldsCk)
        val1  = tryGetField(s, fldsCk{k});
        if isempty(val1), continue; end
        for j = 2:nFile
            % At least some GE ImageOrientationPatient can have diff of 1e-6
            if any(abs(val1 - tryGetField(h{i}{j}, fldsCk{k})) > 1e-4)
                errorLog(['Inconsistent ''' fldsCk{k} ''' for ' series '. Series skipped.']);
                keep(i) = 0;
                break;
            end
        end
        if ~keep(i), break; end % skip
    end
    if ~keep(i), continue; end
    if ~isempty(csa_header(s, 'NumberOfImagesInMosaic')), continue; end

    % this won't catch all missing files, but catch most cases.
    iSL = xform_mat(s); iSL = iSL(3);
    a = zeros(1, nFile);
    for j = 1:nFile, a(j) = h{i}{j}.ImagePositionPatient(iSL); end
    nSL = sum(diff(sort(a)) > 1e-4) + 1; % allow minor diff
    if mod(nFile, nSL)>0
        errorLog(['Missing file(s) detected for ' series '. Series skipped.']);
        keep(i) = 0; continue; % skip
    end
    
    fld = 'LocationsInAcquisition';
    if isfield(s, fld) && s.(fld)~=nSL % warning but override it
        errorLog(sprintf(['%s: the number of slices (%g) in dicom header ' ...
            'seems wrong. nSL = %g is used. Please check the result.'], ...
            series, s.(fld), nSL));
    end
    h{i}{1}.(fld) = uint16(nSL); % best way for nSL?
   
    if strncmpi(s.Manufacturer, 'GE', 2)
        % re-order slices within volume by slice locations: not tested
        [foo, ind] = sort(a(1:nSL)); % slice locations for first volume
        if isequal(ind, 1:nSL) || ~isequal(ind, nSL:-1:1), continue; end
        if foo(1) ~= a(1), ind = ind(nSL:-1:1); end
        nVol = nFile / nSL;
        inc = repmat((0:nVol-1)*nSL, nSL, 1);
        ind = repmat(ind, 1, nVol) + inc(:)';
        h{i} = h{i}(ind); % sorted by slice locations
    elseif strncmpi(s.Manufacturer, 'Philips', 7)
        % re-order as Slice then Volume if Dim3IsVolume
        if abs(a(1)-a(2)) > 1e-4, continue; end % different slices
        nVol = nFile / nSL;
        ind = reshape(1:nFile, [nVol nSL])';
        h{i} = h{i}(ind(:));
        h{i}{1}.Dim3IsVolume = true; % not needed, info only
    end
end

ind = find(isMoCo); % decide MoCo after checking all series
for i = 1:length(ind)
    if MoCo==1 && keep(ind(i)-1) % in case original skipped, keep MOCO
        keep(ind(i)) = 0; continue; % skip MOCO
    elseif MoCo==2 && keep(ind(i)) % in case MOCO skipped, keep original
        keep(ind(i)-1) = 0; % skip previous series (original)
    end
end
h = h(keep); % remove all unwanted series once

%% Generate unique result file names
% Unique names are in format of SeriesDescription_s007. Special cases are: 
%  for phase image, such as field_map phase, append '_phase' to the name;
%  for MoCo series, append '_MoCo' to the name if both series are created.
nRun = length(h); % update it, since we have removed some
if nRun<1
    errorLog(['No valid series found for ' subj]);
    return;
end
rNames = cell(1, nRun);
for i = 1:nRun
    s = h{i}{1};
    a = strtrim(ProtocolName(s));
    a(~isstrprop(a, 'alphanum')) = '_'; % make str valid for field name
    while true % remove repeated underscore
        ind = strfind(a, '__');
        if isempty(ind), break; end
        a(ind) = '';
    end
    if isType(s, '\P\') || strcmpi(tryGetField(s, 'ComplexImageComponent', ''), 'PHASE')
        a = [a '_phase']; % phase image, rely on 1st vol is mag for MIXED series
    end
    if MoCo==0 && isType(s, '\MOCO\'), a = [a '_MoCo']; end
    sN = s.SeriesNumber;
    if sN>100 && strncmp(s.Manufacturer, 'Philips', 7)
        sN = tryGetField(s, 'AcquisitionNumber', floor(sN/100));
    end
    rNames{i} = sprintf('%s_s%03g', a, sN);
end
rNames = genvarname(rNames); % add 'x' if started with a digit, and more

% After following sort, we need to compare only neighboring names. Remove
% _s007 if there is no conflict. Have to ignore letter case for Windows & MAC
fnames = rNames; % copy it, keep letter cases
[rNames, iRuns] = sort(lower(rNames)); 
for i = 1:nRun
    a = rNames{i}(1:end-5); % remove _s003
    % no conflict with both previous and next name
    if nRun==1 || ... % only one run
         (i==1    && ~strcmpi(a, rNames{2}(1:end-5))) || ... % first
         (i==nRun && ~strcmpi(a, rNames{i-1}(1:end-5))) || ... % last
         (i>1 && i<nRun && ~strcmpi(a, rNames{i-1}(1:end-5)) ...
                        && ~strcmpi(a, rNames{i+1}(1:end-5))); % middle ones
        fnames{iRuns(i)}(end+(-4:0)) = [];
    end
end
% fmtStr = sprintf(' %%-%gs %%4g\n', max(cellfun(@length, fnames))+6);
fmtStr = sprintf(' %%-%gs %%gx%%gx%%gx%%g\n', max(cellfun(@length, fnames))+6);

%% Now ready to convert nii series by series
fprintf('Converting %g series into %g-D %s: subject %s\n', ...
            nRun, 4-rst3D, ext, subj);
for i = 1:nRun
    nFile = length(h{i});
    h{i}{1}.isDTI = isDTI(h{i}{1});
    s = h{i}{1};
    img = dicm_img(s, 0); % initialize with proper data type and img size

    if nFile > 1
        h{i}{1}.LastFile = h{i}{nFile}; % store partial last header into 1st
        n = ndims(img);
        if n == 2
            img(:, :, 2:nFile) = 0; % pre-allocate
            for j = 2:nFile, img(:,:,j) = dicm_img(h{i}{j}, 0); end
        elseif n == 3 % SamplesPerPixel>1 is the only case I know for now
            img(:, :, :, 2:nFile) = 0; % pre-allocate
            for j = 2:nFile, img(:,:,:,j) = dicm_img(h{i}{j}, 0); end
        else % err out, likely won't work for other series
            error('dicm2nii can''t deal with multiple %g-dim dicom image', n);
        end
    elseif size(img,3) == 1 % single gray file (multi-frame dicm, par or brik)
        img = permute(img, [1 2 4 3]); % put frames into dim3
    end
    
    nSL = csa_header(s, 'NumberOfImagesInMosaic');
    if ~isempty(nSL) % SIEMENS mosaic
        img = mos2vol(img, nSL); % mosaic to volume
    elseif ndims(img)==4 && tryGetField(s, 'Dim3IsVolume', false)
        img = permute(img, [1 2 4 3]);
    elseif ndims(img) == 3 % may need to reshape to 4D
        nSL = double(tryGetField(s, 'LocationsInAcquisition'));        
        if nSL>1
            dim = size(img);
            dim(3:4) = [nSL dim(3)/nSL]; % verified integer earlier
            if nFile==1 && tryGetField(s, 'Dim3IsVolume', false)
                % for PAR and single multiframe dicom
                img = reshape(img, dim([1 2 4 3]));
                img = permute(img, [1 2 4 3]);
            else
                img = reshape(img, dim);
            end
        end
    end
    % This 'if' block takes care of SamplesPerPixel>1
    if tryGetField(s, 'SamplesPerPixel', 1) > 1
        img = permute(img, [1 2 4:8 3]); % put RGB into dim8 for nii_tool
    end
    dim = size(img);
    if numel(dim)<3, dim(3) = 1; end
    fld = 'NumberOfTemporalPositions';
    if ~isfield(s, fld) && numel(dim)>3 && dim(4)>1, h{i}{1}.(fld) = dim(4); end
        
    % Store GE slice timing. No slice order info for Philips at all!
    if isfield(s, 'RTIA_timer')
        t = zeros(dim(3), 1);
        for j = 1:dim(3), t(j) = tryGetField(h{i}{j}, 'RTIA_timer', nan); end
        if ~all(diff(t)==0), h{i}{1}.SliceTiming = t/10; end % in ms
    end
    
    % Store motion parameters for MoCo series (assume it is mosaic)
    if all(isfield(s, {'RBMoCoTrans' 'RBMoCoRot'}))
        trans = zeros(nFile, 3);
        rotat = zeros(nFile, 3);
        for j = 1:nFile
            trans(j,:) = tryGetField(h{i}{j}, 'RBMoCoTrans', [0 0 0]);
            rotat(j,:) = tryGetField(h{i}{j}, 'RBMoCoRot', [0 0 0]);
        end
        h{i}{1}.RBMoCoTrans = trans;
        h{i}{1}.RBMoCoRot = rotat;
    end
    
    if isa(img, 'uint16') && max(img(:))<32768
        img = int16(img); % use int16 if lossless. seems always true
    end
    
    nii = nii_tool('init', img);
    fname = [dataFolder fnames{i}];

    % get bval & bvec in dicom image reference
    if s.isDTI, [h{i}, nii] = get_dti_para(h{i}, nii); end
    
    [nii, h{i}{1}] = set_nii_header(nii, h{i}{1}); % set most nii header

    % Save FSL bval and bvec files after computing bve in set_nii_header
    if s.isDTI, save_dti_para(h{i}{1}, fname); end
    
    [nii, niiP] = split_philips_phase(nii, s);
    if ~isempty(niiP)
        fprintf(fmtStr, [fnames{i} '_phase'], niiP.hdr.dim(2:5));
        nii_tool('save', niiP, [fname '_phase' ext], rst3D); % save phase nii
    end
    
    fprintf(fmtStr, fnames{i}, nii.hdr.dim(2:5)); % show info and progress
    nii_tool('save', nii, [fname ext], rst3D);
    h{i} = h{i}{1}; % keep 1st dicm header only
end

h = cell2struct(h, fnames, 2); % convert into struct
fname = [dataFolder 'dcmHeaders.mat'];
if exist(fname, 'file') % if file exists, we update fields only
    S = load(fname);
    for i = 1:length(fnames), S.h.(fnames{i}) = h.(fnames{i}); end
    h = S.h; %#ok
end
save(fname, 'h', '-v7'); % -v7 better compatibility
fprintf('Elapsed time by dicm2nii is %.1f seconds\n\n', toc);
if ~isempty(unzip_cmd), rmdir(dcmFolder, 's'); end % delete tmp dicom folder
return;

%% Subfunction: return folder name for a file name
function folder = folderFromFile(fname)
folder = fileparts(fname);
if isempty(folder), folder = pwd; end

%% Subfunction: return SeriesDescription
function name = ProtocolName(s)
name = tryGetField(s, 'SeriesDescription');
if isempty(name), name = tryGetField(s, 'ProtocolName', ''); end
if strncmpi(s.Manufacturer, 'SIEMENS', 7) || isempty(name)
    name = tryGetField(s, 'ProtocolName', name); 
end
if isempty(name), [~, name] = fileparts(s.Filename); end

%% Subfunction: return true if any of keywords is in s.ImageType
function tf = isType(s, keywords)
if ischar(keywords) % single keyword
    tf = ~isempty(strfind(s.ImageType, keywords));
    return;
end
for i = 1:length(keywords)
    tf = ~isempty(strfind(s.ImageType, keywords{i}));
    if tf, return; end
end

%% Subfunction: return true if series is DTI
function tf = isDTI(s)
tf = isType(s, '\DIFFUSION'); % Siemens, Philips
if tf, return; end
if strncmp(s.Manufacturer, 'GE', 2) % not labeled as /DIFFISION
    tf = tryGetField(s, 'DiffusionDirection', 0)>0;
elseif strncmpi(s.Manufacturer, 'Philips', 7)
    tf = strcmp(tryGetField(s, 'MRSeriesDiffusion', 'N'), 'Y');
else % Some Siemens DTI are not labeled as \DIFFUSION
    tf = ~isempty(csa_header(s, 'B_value'));
end
        
%% Subfunction: get field if exist, return default value otherwise
function val = tryGetField(s, field, dftVal)
if isfield(s, field), val = s.(field); 
elseif nargin>2, val = dftVal;
else val = [];
end

%% Subfunction: Set most nii header and re-orient img
function [nii, s] = set_nii_header(nii, s)
TR = tryGetField(s, 'RepetitionTime'); % in ms
if isempty(TR), TR = tryGetField(s, 'TemporalResolution', 2000); end
dim = nii.hdr.dim(2:4); % image dim, set by nii_tool according to img

% Transformation matrix: most important feature for nii
[ixyz, R, pixdim, revSL] = xform_mat(s, dim);
R(1:2,:) = -R(1:2,:); % dicom LPS to nifti RAS, xform matrix before reorient

% dim_info byte: freq_dim, phase_dim, slice_dim low to high, each 2 bits
[phPos, fps_bits] = phaseDirection(s); % phPos relative to image in FSL feat!

% sign change at slice direction due to possibly reversed slice order
if revSL && isfield(s, 'bvec'), s.bvec(:,3) = -s.bvec(:,3); end

% Reorient if MRAcquisitionType==3D || isDTI && dim(3)>1
% If FSL etc can read dim_info for STC, we can always reorient.
[~, perm] = sort(ixyz); % may permute 3 dimensions in this order
if (strcmp(tryGetField(s, 'MRAcquisitionType'), '3D') || s.isDTI) && ...
        dim(3)>1 && (~isequal(perm, 1:3)) % skip if already standard view
    R(:, 1:3) = R(:, perm); % xform matrix after perm
    fps_bits = fps_bits(perm);
    ixyz = ixyz(perm); % 1:3 after re-orient
    dim = dim(perm);
    pixdim = pixdim(perm);
    nii.hdr.dim(2:4) = dim;
    nii.img = permute(nii.img, [perm 4:8]); % permute img after flip
    if isfield(s, 'bvec'), s.bvec = s.bvec(:,perm); end
end
iPhase = find(fps_bits==4); % now it is axis index for phase_dim

% Flip image to make first axis negative and other two positive
ind4 = ixyz + [0 4 8]; % index in 4xN matrix
flip = R(ind4)<0; % flip an axis if true
flip(1) = ~flip(1); % first axis negative: comment this to make all positive
rotM = diag([sign(0.5-flip) 1]); % 1 or -1 on diagnal
rotM(1:3, 4) = (dim-1) .* flip;
R = R / rotM; % xform matrix after flip
for k = 1:3, if flip(k), nii.img = flipdim(nii.img, k); end; end
if flip(iPhase), phPos = ~phPos; end
if isfield(s, 'bvec'), s.bvec(:, flip) = -s.bvec(:, flip); end

frmCode = tryGetField(s, 'TemplateSpace', 1); % 1: SCANNER_ANAT
nii.hdr.sform_code = frmCode;
nii.hdr.srow_x = R(1,:);
nii.hdr.srow_y = R(2,:);
nii.hdr.srow_z = R(3,:);

nii.hdr.qform_code = frmCode;
nii.hdr.qoffset_x = R(1,4);
nii.hdr.qoffset_y = R(2,4);
nii.hdr.qoffset_z = R(3,4);

R = R(1:3, 1:3); % for quaternion
R = R ./ repmat(sqrt(sum(R.^2)),3,1); % normalize
[q, nii.hdr.pixdim(1)] = dcm2quat(R); % 3x3 dir cos matrix to quaternion
nii.hdr.quatern_b = q(2);
nii.hdr.quatern_c = q(3);
nii.hdr.quatern_d = q(4);

str = tryGetField(s, 'ImageComments');
if isType(s, '\MOCO\'), str = ''; end % useless for MoCo
foo = tryGetField(s, 'StudyComments');
if ~isempty(foo), str = [str ';' foo]; end
str = [str ';' strtok(s.Manufacturer)];
foo = tryGetField(s, 'ProtocolName');
if ~isempty(foo), str = [str ';' foo]; end
nii.hdr.aux_file = str; % char[24], info only
seq = asc_header(s, 'tSequenceFileName'); % like '%SiemensSeq%\ep2d_bold'
[~, seq] = strtok(seq, '\'); seq = strtok(seq, '\'); % like 'ep2d_bold'
if isempty(seq), seq = tryGetField(s, 'ScanningSequence'); end
id = tryGetField(s, 'PatientID');
nii.hdr.db_name = id; % char[18], optional
nii.hdr.intent_name = seq; % char[16], meaning of the data

% save some useful info in descrip: info only
foo = tryGetField(s, 'AcquisitionDateTime');
if isempty(foo) 
    foo = tryGetField(s, 'AcquisitionDate');
    foo = [foo tryGetField(s, 'AcquisitionTime')];
end
descrip = sprintf('time=%s;', foo(1:min(18,end))); 
TE0 = asc_header(s, 'alTE[0]')/1000; % s.EchoTime stores only 1 TE
dTE = abs(asc_header(s, 'alTE[1]')/1000 - TE0); % TE difference
if isempty(TE0), TE0 = tryGetField(s, 'EchoTime'); end % GE, philips
if isempty(dTE) && tryGetField(s, 'NumberOfEchoes', 1)>1
    dTE = tryGetField(s, 'SecondEchoTime') - TE0; % need to update
end
if ~isempty(dTE)
    descrip = sprintf('dTE=%.4g;%s', dTE, descrip);
    s.deltaTE = dTE;
elseif ~isempty(TE0)
    descrip = sprintf('TE=%.4g;%s', TE0, descrip);
end

% Get dwell time, slice timing info
if ~strcmp(tryGetField(s, 'MRAcquisitionType'), '3D')
    hz = csa_header(s, 'BandwidthPerPixelPhaseEncode');
    dwell = 1000 ./ hz / dim(iPhase); % in ms
    if isempty(dwell) % true for syngo MR 2004A
        % ppf = [1 2 4 8 16] represent [4 5 6 7 8] 8ths PartialFourier
        % ppf = asc_header(s, 'sKSpace.ucPhasePartialFourier');
        lns = asc_header(s, 'sKSpace.lPhaseEncodingLines');
        dur = csa_header(s, 'SliceMeasurementDuration');
        dwell = dur ./ lns; % ./ (log2(ppf)+4) * 8;
    end
    if isempty(dwell) % next is not accurate, so as last resort
        dur = csa_header(s, 'RealDwellTime') * 1e-6; % ns to ms
        dwell = dur * asc_header(s, 'sKSpace.lBaseResolution');
    end
    if isempty(dwell)
        dwell = double(tryGetField(s, 'EffectiveEchoSpacing')) / 1000; % GE
    end
    % http://www.spinozacentre.nl/wiki/index.php/NeuroWiki:Current_developments
    if isempty(dwell) % Philips
        wfs = tryGetField(s, 'WaterFatShift');
        epiFactor = tryGetField(s, 'EPIFactor');
        dwell = wfs ./ (434.215 * (double(epiFactor)+1)) * 1000;
    end
    if ~isempty(dwell)
        if s.isDTI
            readout = dwell * dim(iPhase) / 1000; % in sec
            descrip = sprintf('readout=%.3g;%s', readout, descrip);
            s.ReadoutSeconds = readout;
        else
            descrip = sprintf('dwell=%.3g;%s', dwell, descrip);
            s.EffectiveEPIEchoSpacing = dwell;
        end
    end
    
    t = csa_header(s, 'MosaicRefAcqTimes'); % in ms
    % MosaicRefAcqTimes for first vol may be wrong for Siemens MB
    if ~isempty(t) && isfield(s, 'LastFile') % pity: no MB flag in dicom
        dict = dicm_dict(s.Manufacturer, 'MosaicRefAcqTimes');
        s2 = dicm_hdr(s.LastFile.Filename, dict);
        t = s2.MosaicRefAcqTimes; % to be safe, use the last file
    end
    if isempty(t), t = tryGetField(s, 'SliceTiming'); end % for GE
    
    if numel(t)>1 % 1+ slices
        t1 = sort(t);
        dur = mean(diff(t1));
        dif = mean(diff(t));
        if dur==0 || (t1(end)>TR), sc = 0; % no useful info, or bad timing MB
        elseif t1(1) == t1(2), sc = 7; % timing available MB, made-up code 7
        elseif abs(dif-dur)<1e-3, sc = 1; % ascending
        elseif abs(dif+dur)<1e-3, sc = 2; % descending
        elseif t(1)<t(3) % ascending interleaved
            if t(1)<t(2), sc = 3; % odd slices first
            else sc = 5; % Siemens even number of slices
            end
        elseif t(1)>t(3) % descending interleaved
            if t(1)>t(2), sc = 4;
            else sc = 6; % Siemens even number of slices
            end
        else sc = 0; % unlikely to reach
        end
 
        if flip(fps_bits == 16) % slices flipped
            if sc>0, sc = sc+mod(sc,2)*2-1; end % 1<->2, 3<->4, 5<->6
            t = t(end:-1:1);
        end
        t = t - min(t); % it may be relative to 1st slice
        if sc>0, s.SliceTiming = 0.5 - t/TR; end % as for FSL custom timing
        nii.hdr.slice_code = sc;
    end
    nii.hdr.slice_end = dim(3)-1; % 0-based, slice_start default to 0

    dur = min(diff(sort(t))); % 2.5ms error Siemens, dur = 0 for MB
    if isempty(dur) % in case MosaicRefAcqTimes is not available
        delay = asc_header(s, 'lDelayTimeInTR')/1000; % in ms now
        if isempty(delay), delay = 0; end
        dur = (TR-delay)/dim(3); 
    end
    nii.hdr.slice_duration = dur/1000; 
end

% data slope and intercept: apply to img if no rounding error 
nii.hdr.scl_slope = 1; % default scl_inter is 0
if any(isfield(s, {'RescaleSlope' 'RescaleIntercept'}))
    slope = tryGetField(s, 'RescaleSlope', 1); 
    inter = tryGetField(s, 'RescaleIntercept', 0); 
    val = sort([max(nii.img(:)) min(nii.img(:))] * slope + inter);
    dClass = class(nii.img);
    if isa(nii.img, 'float') || (mod(slope,1)==0 && mod(inter,1)==0 ... 
            && val(1)>=intmin(dClass) && val(2)<=intmax(dClass))
        nii.img = nii.img * slope + inter; % apply to img if no rounding
    else
        nii.hdr.scl_slope = slope;
        nii.hdr.scl_inter = inter;
    end
end

if isempty(phPos), pm = '?'; elseif phPos, pm = ''; else pm = '-'; end
axes = 'xyz';
phDir = [pm axes(iPhase)];
s.UnwarpDirection = phDir;
descrip = sprintf('phase=%s;%s', phDir, descrip);
s.NiftiCreator = ['dicm2nii.m 20' reviseDate];

nii.hdr.descrip = descrip; % char[80], drop from end if exceed
nii.hdr.dim_info = (1:3) * fps_bits'; % useful for EPI only
nii.hdr.xyzt_units = 10; % mm (2) + seconds (8)
nii.hdr.pixdim(2:5) = [pixdim TR/1000]; % voxSize and TR

% Possible patient position: HFS/HFP/FFS/FFP / HFDR/HFDL/FFDR/FFDL
% Seems dicom takes care of this, and maybe nothing needs to do here.
% patientPos = tryGetField(s, 'PatientPosition', '');

% Save 1st dicom hdr into nii extension
persistent saveExt40
if isempty(saveExt40)
    % FSL 5.0.5 and 5.0.6 have problem with ecode>30. For now, we save ext only
    % if we detect FSL version, and it is not 5.0.5 or 5.0.6
    fid = fopen([getenv('FSLDIR') '/etc/fslversion']);
    if fid < 1
        fslver = '';
    else
        fslver = strtrim(fread(fid, '*char')');
        fclose(fid);
    end
    saveExt40 = ~(isempty(fslver) || ~isempty(strfind(fslver, '5.0.5')) || ...
        ~isempty(strfind(fslver, '5.0.6')));
end

if saveExt40
    fname = [tempname '.mat'];
    save(fname, '-struct', 's', '-v7'); % field as variable
    fid = fopen(fname);
    b = fread(fid, inf, '*uint8'); % data bytes
    fclose(fid);
    delete(fname);
    
    % first 4 bytes (int32) encode real data length. Don't care endian since it
    % will be typecast'ed into int32, rather than fread as int32.
    nii.ext.edata = [typecast(int32(numel(b)), 'uint8')'; b];
    nii.ext.ecode = 40; % Matlab
end

% n = 1; if isfield(nii, 'ext'), n = length(nii.ext)+1; end
% nii.ext(n).ecode = 2; % dicom
% fid = fopen(s.Filename);
% nii.ext(n).edata = fread(fid, s.PixelData.Start, '*uint8');
% fclose(fid);

%% Subfunction, reshape mosaic into volume, remove padded zeros
function vol = mos2vol(mos, nSL)
nMos = ceil(sqrt(nSL)); % always nMos x nMos tiles
[nr, nc, nv] = size(mos); % number of row, col and vol in mosaic

nr = nr / nMos; nc = nc / nMos; % number of row and col in slice
vol = zeros([nr nc nSL nv], class(mos));
for i = 1:nSL
    r =    mod(i-1, nMos) * nr + (1:nr); % 2nd slice is tile(2,1)
    c = floor((i-1)/nMos) * nc + (1:nc);
    % r = floor((i-1)/nMos) * nr + (1:nr); % 2nd slice is tile(1,2)
    % c =    mod(i-1, nMos) * nc + (1:nc);
    vol(:, :, i, :) = mos(r, c, :);
end

%% subfunction: extract bval & bvec, store in 1st header
function [h, nii] = get_dti_para(h, nii)
nSL = nii.hdr.dim(4);
nDir = nii.hdr.dim(5);
bval = nan(nDir, 1);
bvec = nan(nDir, 3);
s = h{1};

if tryGetField(s, 'IsPhilipsPAR', false) || tryGetField(s, 'IsAFNIHEAD', false)
    bval = s.B_value;
    bvec = tryGetField(s, 'DiffusionGradientDirection', nan(nDir, 3));
elseif isfield(s, 'PerFrameFunctionalGroupsSequence')
    fld = 'PerFrameFunctionalGroupsSequence';
    if tryGetField(s, 'Dim3IsVolume', false), iDir = 1:nDir;
    else iDir = 1:nSL:nSL*nDir;
    end
    dict = dicm_dict(s.Manufacturer, {fld 'B_value' 'MRDiffusionSequence' ...
        'DiffusionGradientDirectionSequence' 'DiffusionGradientDirection'});
    s2 = dicm_hdr(s.Filename, dict, iDir); % re-read needed frames
    sq = s2.(fld);
    for j = 1:nDir
        item = sprintf('Item_%g', iDir(j));
        try
            a = sq.(item).MRDiffusionSequence.Item_1;
            bval(j) = a.B_value;
            a = a.DiffusionGradientDirectionSequence.Item_1;
            bvec(j,:) = a.DiffusionGradientDirection;
        end
    end
else % multiple files: order already in slices then volumes
    dict = dicm_dict(s.Manufacturer, {'B_value' 'SlopInt_6_9' ...
       'DiffusionDirectionX' 'DiffusionDirectionY' 'DiffusionDirectionZ'});
    iDir = (0:nDir-1) * length(h)/nDir + 1; % could be mosaic 
    for j = 1:nDir % no these values for 1st file of each excitation
        s2 = h{iDir(j)};
        vec = tryGetField(s2, 'DiffusionGradientDirection');
        if isempty(vec) % GE, Philips
            s2 = dicm_hdr(s2.Filename, dict);
            vec(1) = tryGetField(s2, 'DiffusionDirectionX', 0);
            vec(2) = tryGetField(s2, 'DiffusionDirectionY', 0);
            vec(3) = tryGetField(s2, 'DiffusionDirectionZ', 0);
        end
        bvec(j,:) = vec;
        
        val = tryGetField(s2, 'B_value');
        if isempty(val) && isfield(s2, 'SlopInt_6_9') % GE
            val = s2.SlopInt_6_9(1);
        end
        if isempty(val), val = 0; end % B0
        bval(j) = val;
    end
end

if all(isnan(bval)) && all(isnan(bvec(:)))
    errorLog(['Failed to get DTI parameters: ' ProtocolName(s)]);
    return; 
end
bval(isnan(bval)) = 0;
bvec(isnan(bvec)) = 0;

if strncmpi(s.Manufacturer, 'Philips', 7)
    if max(sum(bvec.^2, 2)) > 2 % guess in degree
        for j = 1:nDir, bvec(j,:) = ang2vec(bvec(j,:)); end % deg to dir cos mat
        errorLog(['Please validate bvec (direction in deg): ' ProtocolName(s)]);
    end
    if bval(nDir)~=0 && all(abs(bvec(nDir,:))<1e-4)
        % Remove last vol if it is computed ADC
        bval(nDir) = [];
        bvec(nDir,:) = [];
        nii.img(:,:,:,nDir) = [];
        nDir = nDir - 1;
        nii.hdr.dim(5) = nDir;
    end
end

h{1}.B_value = bval; % store all into header of 1st file
h{1}.DiffusionGradientDirection = bvec; % original from dicom

% http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:DICOM_for_DWI_and_DTI
if strncmp(s.Manufacturer, 'GE', 2) % GE bvec already in image reference
    % Following sign change is based on FSL result. non-axial slice not tested
    if strcmp(tryGetField(s, 'InPlanePhaseEncodingDirection'), 'ROW')
        bvec = bvec(:, [2 1 3]);
        bvec(:, 2:3) = -bvec(:, 2:3);
    else
        bvec(:, 1:2) = -bvec(:, 1:2);
    end
else % Siemens/Philips
    R = reshape(s.ImageOrientationPatient, 3, 2);
    R(:,3) = null(R'); % det(R) is always 1
    bvec = bvec * R; % dicom plane to image plane
end

h{1}.bvec = bvec; % computed bvec

%% subfunction: save bval & bvec files
function save_dti_para(s, fname)
if isfield(s, 'B_value')
    fid = fopen([fname '.bval'], 'w');
    fprintf(fid, '%g ', s.B_value); % one row
    fclose(fid);
end

if isfield(s, 'bvec')
    str = repmat('%.6g ', 1, size(s.bvec,1));
    fid = fopen([fname '.bvec'], 'w');
    fprintf(fid, [str '\n'], s.bvec); % 3 rows by # direction cols
    fclose(fid);
end

%% Subfunction: convert rotation angles to vector
function vec = ang2vec(ang)
% do the same as in philips_par: not sure it is right
ca = cosd(ang); sa = sind(ang);
rx = [1 0 0; 0 ca(1) -sa(1); 0 sa(1) ca(1)]; % standard 3D rotation
ry = [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)];
rz = [ca(3) -sa(3) 0; sa(3) ca(3) 0; 0 0 1];
R = rx * ry * rz;
% [~, iSL] = max(abs(R(:,3)));
% if iSL == 1 % Sag
%     R(:,[1 3]) = -R(:,[1 3]);
%     R = R(:, [2 3 1]);
% elseif iSL == 2 % Cor
%     R(:,3) = -R(:,3);
%     R = R(:, [1 3 2]);
% end
% http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
vec = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
vec = vec / sqrt(sum(vec.^2));

%% Subfunction, return a parameter from CSA header
function val = csa_header(s, key)
try val = s.CSAImageHeaderInfo.(key); 
catch, val = [];
end

%% Subfunction, Convert 3x3 direction cosine matrix to quaternion
% Simplied from Quaternions by Przemyslaw Baranski 
function [q, proper] = dcm2quat(R)
proper = sign(det(R)); % always -1 if flip_dim1
if proper<0, R(:,3) = -R(:,3); end

q = zeros(1, 4);
q(1) = sqrt(1 + R(1,1) + R(2,2) + R(3,3)) / 2;
q(2) = sqrt(1 + R(1,1) - R(2,2) - R(3,3)) / 2;
q(3) = sqrt(1 - R(1,1) + R(2,2) - R(3,3)) / 2;
q(4) = sqrt(1 - R(1,1) - R(2,2) + R(3,3)) / 2;
if ~isreal(q(1)), q(1) = 0; end % if trace(R)+1<0, zero it
[m, ind] = max(q);

switch ind
    case 1
        q(2) = (R(3,2) - R(2,3)) /m/4;
        q(3) = (R(1,3) - R(3,1)) /m/4;
        q(4) = (R(2,1) - R(1,2)) /m/4;
    case 2
        q(1) = (R(3,2) - R(2,3)) /m/4;
        q(3) = (R(1,2) + R(2,1)) /m/4;
        q(4) = (R(3,1) + R(1,3)) /m/4;
    case 3
        q(1) = (R(1,3) - R(3,1)) /m/4;
        q(2) = (R(1,2) + R(2,1)) /m/4;
        q(4) = (R(2,3) + R(3,2)) /m/4;
    case 4
        q(1) = (R(2,1) - R(1,2)) /m/4;
        q(2) = (R(3,1) + R(1,3)) /m/4;
        q(3) = (R(2,3) + R(3,2)) /m/4;
end
if q(1)<0, q = -q; end % as MRICron

%% Subfunction: get dicom xform matrix and related info
function [ixyz, R, pixdim, flip] = xform_mat(s, dim)
R = reshape(s.ImageOrientationPatient, 3, 2);
R(:, 3) = null(R');
foo = abs(R);
[~, ixyz] = max(foo); % orientation info: perm of 1:3
if ixyz(2) == ixyz(1), foo(ixyz(2),2) = 0; [~, ixyz(2)] = max(foo(:,2)); end
if ixyz(3) == ixyz(1), foo(ixyz(3),3) = 0; [~, ixyz(3)] = max(foo(:,3)); end
if ixyz(3) == ixyz(2), foo(ixyz(3),3) = 0; [~, ixyz(3)] = max(foo(:,3)); end
if nargout<2, return; end
flip = false; 

iSL = ixyz(3); % 1/2/3 for Sag/Cor/Tra slice
if isfield(s, 'LastFile') && (dim(1) == s.Columns) % most reliable way
    thk = s.LastFile.ImagePositionPatient - s.ImagePositionPatient;
    thk = abs(thk(iSL) / R(iSL,3)) / (dim(3)-1);
else % mainly mosaic
    thk = tryGetField(s, 'SpacingBetweenSlices');
    if isempty(thk), thk = tryGetField(s, 'SliceThickness'); end
    if isempty(thk)
        errorLog(['No slice thickness information found: ' ProtocolName(s)]);
        thk = s.PixelSpacing(1); % guess it is cubic
    end
end

pixdim = [s.PixelSpacing(:)' thk];
R = R * diag(pixdim); % apply vox size
% Next is almost dicom xform matrix, except mosaic trans and unsure slice_dir
R = [R s.ImagePositionPatient; 0 0 0 1];

% rest are former: R = verify_slice_dir(R, s, dim, iSL)
if dim(3)<2, return; end % don't care direction for single slice
pos = []; % SliceLocation for last or center slice we try to retrieve

if s.Columns > dim(1) % Siemens mosaic: use Columns since no transpose to img
    R(:,4) = R * [ceil(sqrt(dim(3))-1)*dim(1:2)/2 0 1]'; % transform back
    vec = csa_header(s, 'SliceNormalVector'); % mosaic has this
    if ~isempty(vec) % 100% exist for all tested data
        flip = sign(vec(iSL)) ~= sign(R(iSL,3));
        if flip, R(:,3) = -R(:,3); end
        return;
    end
    
    % Use 2nd slice so no worry for revNumbering.
    % This wont solve nSL=2 case, but here is almost unreachable anyway. 
    if isfield(s, 'CSASeriesHeaderInfo') && dim(3)>2
        % sSliceArray.asSlice[0].sPosition.dSag/Cor/Tra
        ori = ['Sag'; 'Cor'; 'Tra']; % 1/2/3
        pos = asc_header(s, ['sSliceArray.asSlice[1].sPosition.d' ori(iSL,:)]);
        if ~isempty(pos)
            pos = [R(iSL, 1:2) pos] * [-dim(1:2)/2 1]'; % 2nd slice location
        end
    end
end

% s.LastFile works for most GE, Philips and all non-mosaic Siemens data
if isempty(pos) && isfield(s, 'LastFile')
    pos = tryGetField(s.LastFile, 'ImagePositionPatient');
    if ~isempty(pos)
        pos = pos(iSL);
        if pos == s.ImagePositionPatient(iSL), pos = []; end % ignore same slice
    end
end

% May be useful for Philips dicom: use volume centre info
if isempty(pos) && isfield(s, 'Stack')
    ori = ['RL'; 'AP'; 'FH']; % x y z
    pos = tryGetField(s.Stack.Item_1, ['MRStackOffcentre' ori(iSL,:)]);
    if ~isempty(pos)
        pos = [R(iSL, 1:2) pos] * [-dim(1:2)/2 1]'; % mid-slice location
    end
end

% GE: LastScanLoc is always different from the slice in 1st file
if isempty(pos) && isfield(s, 'LastScanLoc')
    pos = s.LastScanLoc;
    if iSL<3, pos = -pos; end % LastScanLoc uses RAS convention!
end

pos1 = R(iSL, 3:4) * [dim(3)-1 1]'; % last SliceLocation based on current R
if ~isempty(pos) % have SliceLocation for last/center slice
    flip = (pos>R(iSL,4)) ~= (pos1>R(iSL,4)); % same direction?
else % we do some guess work and warn user
    errorLog(['Please check whether slices are flipped: ' ProtocolName(s)]);
    pos2 = R(iSL, 3:4) * [1-dim(3) 1]'; % opposite slice direction
    % if pos1 is larger than the other dir, and is way outside head
    flip = all(abs(pos1) > [abs(pos2) 150]); % arbituary 150 mm
end
if flip, R(:,3) = -R(:,3); end % change to opposite direction

%% Subfunction: get a parameter in CSA series ASC header.
function val = asc_header(s, key)
val = []; 
fld = 'CSASeriesHeaderInfo';
if ~isfield(s, fld), return; end
if isfield(s.(fld), key)
    val = s.(fld).(key);
    return;
end
if isfield(s.(fld), 'MrPhoenixProtocol')
    str = s.(fld).MrPhoenixProtocol;
elseif isfield(s.(fld), 'MrProtocol') % older version dicom
    str = s.(fld).MrProtocol;
else
    str = char(s.(fld)');
    k0 = strfind(str, '### ASCCONV BEGIN ###');
    k  = strfind(str, '### ASCCONV END ###');
    str = str(k0:k); % avoid key before BEGIN and after END
end
k = strfind(str, [char(10) key]); % start with new line: safer
if isempty(k), return; end
str = strtok(str(k(1):end), char(10)); % the line
[~, str] = strtok(str, '='); % '=' and the vaule
str = strtrim(strtok(str, '=')); % remvoe '=' and space 

if strncmp(str, '""', 2) % str parameter
    val = str(3:end-2);
elseif strncmp(str, '"', 1) % str parameter for version like 2004A
    val = str(2:end-1);
elseif strncmp(str, '0x', 2) % hex parameter, convert to decimal
    val = sscanf(str(3:end), '%x', 1);
else % decimal
    val = str2double(str);
end

%% Subfunction: return matlab decompress command if the file is compressed
function func = compress_func(fname)
func = '';
fid = fopen(fname);
if fid<0, return; end
sig = fread(fid, 4, '*uint8')';
fclose(fid);
if isequal(sig, [80 75 3 4]) % zip file
    func = 'unzip';
elseif isequal(sig, [31 139 8 0]) % tgz, tar
    func = 'untar';
end
% ! "c:\program Files (x86)\7-Zip\7z.exe" x -y -oF:\tmp\ F:\zip\3047ZL.zip

%% Subfuction: for GUI subfunctions
function dicm2nii_gui(cmd, hs)
if nargin<2, hs = guidata(gcbo); end
drawnow;
switch cmd
    case 'do_convert'
        if get(hs.srcType, 'Value') > 2 % dicom, PAR, HEAD files
            src = get(hs.src, 'UserData');
        else
            src = get(hs.src, 'String');
        end
        dst = get(hs.dst, 'String');
        if isempty(src) || isempty(dst)
            str = 'Dicom source and Result folder must be specified';
            errordlg(str, 'Error Dialog');
            return;
        end
        rstFmt = (get(hs.rstFmt, 'Value') - 1) * 2; % 0 or 2
        if get(hs.gzip,  'Value'), rstFmt = rstFmt + 1; end % 1 or 3 
        if get(hs.rst3D, 'Value'), rstFmt = rstFmt + 4; end % 4 to 7
        mocoOpt = get(hs.mocoOpt, 'Value') - 1;
        subjName = strtrim(get(hs.subjName, 'String'));
        if get(hs.subjPop, 'Value')==1, subjName = ''; end
        set(hs.convert, 'Enable', 'off', 'string', 'Conversion in progress');
        drawnow;
        
        para.srcType = get(hs.srcType, 'Value');
        para.rstFmt = get(hs.rstFmt, 'Value');
        para.rst3D = get(hs.rst3D, 'Value');
        para.gzip = get(hs.gzip, 'Value');
        para.mocoOpt = get(hs.mocoOpt, 'Value');
        para.subjPop = get(hs.subjPop, 'Value');
        para.src = get(hs.src, 'String');
        para.dst = get(hs.dst, 'String');
        fname = [fileparts(which(mfilename)) '/dicm2nii_gui_para.mat'];
        try save(fname, 'para'); end
        
        try
            dicm2nii(src, dst, rstFmt, mocoOpt, subjName);
        catch me
            try set(hs.convert, 'Enable', 'on', 'String', 'Start conversion'); end
            commandwindow;
            rethrow(me);
        end
        try %#ok<*TRYNC>
            set(hs.convert, 'Enable', 'on', 'String', 'Start conversion');
        end
    case 'srcType'
        i = get(hs.srcType, 'Value');
        txt = {'Dicom folder' 'Zip / tar file' 'Dicom files' 'PAR files' 'HEAD files'};
        set(hs.srcTxt, 'String' , txt{i});
        srcdir = isdir(get(hs.src, 'String'));
        if (i==1 && ~srcdir) || (i>1 && srcdir)
            set(hs.src, 'String', '');
        end
    case 'dstDialog'
        folder = get(hs.dst, 'String'); % current folder
        if ~isdir(folder), folder = get(hs.src, 'String'); end
        if ~isdir(folder), folder = fileparts(folder); end
        if ~isdir(folder), folder = pwd; end
        dst = uigetdir(folder, 'Select a folder to save data files');
        if isnumeric(dst), return; end
        set(hs.dst, 'String' , dst);
    case 'srcDialog'
        folder = get(hs.src, 'String'); % initial folder
        if ~isdir(folder), folder = fileparts(folder); end
        if ~isdir(folder), folder = pwd; end
        i = get(hs.srcType, 'Value');
        if i == 1 % folder
            src = uigetdir(folder, 'Select a folder containing DICOM files');
            if isnumeric(src), return; end
            set(hs.src, 'UserData', src);
        elseif i == 2 % zip/tgz file
            [src, folder] = uigetfile([folder '/*.zip;*.tgz;*.tar;*.tar.gz'], ...
                'Select the compressed file containing data files');
            if isnumeric(src), return; end
            src = fullfile(folder, src);
            set(hs.src, 'UserData', src);
        elseif i == 3 % dicom files
            [src, folder] = uigetfile([folder '/*.dcm'], ...
                'Select one or more DICOM files', 'MultiSelect', 'on');
            if isnumeric(src), return; end
            src = cellstr(src); % in case only 1 file selected
            src = strcat(folder, filesep, src);
            set(hs.src, 'UserData', src);
            src = src{1};
        elseif i == 4 % PAR/REC
            [src, folder] = uigetfile([folder '/*.PAR'], ...
                'Select one or more PAR files', 'MultiSelect', 'on');
            if isnumeric(src), return; end
            src = cellstr(src); % in case only 1 file selected
            src = strcat(folder, src);
            set(hs.src, 'UserData', src);
            src = src{1};
        elseif i == 5 % HEAD/BRIK
            [src, folder] = uigetfile([folder '/*.HEAD'], ...
                'Select one or more HEAD files', 'MultiSelect', 'on');
            if isnumeric(src), return; end
            src = cellstr(src); % in case only 1 file selected
            src = strcat(folder, src);
            set(hs.src, 'UserData', src);
            src = src{1};
        end
        set(hs.src, 'String' , src);
        dicm2nii_gui('subjPop', hs);
    case 'set_src'
        str = get(hs.src, 'String');
        if ~exist(str, 'file')
            val = dir(str);
            folder = fileparts(str);
            if isempty(val)
                val = get(hs.src, 'UserData');
                if iscellstr(val), val = val{1}; end
                set(hs.src, 'String', val);
                errordlg('Invalid input', 'Error Dialog');
                return;
            end
            str = {val.name};
            str = strcat(folder, filesep, str);
        end
        set(hs.src, 'UserData', str);
        dicm2nii_gui('subjPop', hs);
    case 'set_dst'
        str = get(hs.dst, 'String');
        if isempty(str), return; end
        if ~exist(str, 'file') && ~mkdir(str)
            set(hs.dst, 'String', '');
            errordlg(['Invalid folder name ''' str ''''], 'Error Dialog');
            return;
        end
    case 'subjPop'
        if get(hs.subjPop, 'Value')==2 && ...
                ~isempty(strtrim(get(hs.subjName, 'String')))
            return;
        end
        src = get(hs.src, 'UserData');
        if isempty(src), return; end
        if iscellstr(src)
            fname = src;
        elseif isdir(src)
            fname = dir(src);
            fname([fname.isdir]) = [];
            fname = strcat(src, filesep, {fname.name});
        else
            fname = cellstr(src);
        end
        dict = dicm_dict('', {'PatientName' 'PatientID'});
        for i = 1:length(fname)
            s = dicm_hdr(fname{i}, dict);
            if isempty(s), continue; end
            subj = tryGetField(s, 'PatientName');
            if isempty(subj)
                subj = tryGetField(s, 'PatientID', 'unknown'); 
            end
            set(hs.subjName, 'String', subj);
            break;
        end
    case 'SPMStyle' % turn off compression
        if get(hs.rst3D, 'Value'), set(hs.gzip, 'Value', 0); end
    otherwise
        create_gui;
end

%% Subfuction: create GUI or bring it to front if exists
function create_gui
fh = figure(typecast(uint8('dicm'), 'uint32')); % arbitury integer
if strcmp('dicm2nii_fig', get(fh, 'Tag')), return; end
scrSz = get(0, 'ScreenSize');
set(fh, 'Toolbar', 'none', 'Menubar', 'none', 'Resize', 'off', ...
    'Tag', 'dicm2nii_fig', 'Position', [200 scrSz(4)-500 420 300], ...
    'Name', 'DICOM to NIfTI Converter', 'NumberTitle', 'off');
clr = get(fh, 'color');

str = 'Choose what kind of data source you are using';
uicontrol('Style', 'text', 'Position', [10 250 90 30], ...
    'FontSize', 9, 'HorizontalAlignment', 'left', ...
    'String', 'Source type', 'Background', clr, 'TooltipString', str);
uicontrol('Style', 'popup', 'Background', 'white', 'Tag', 'srcType', ...
    'String', [' Folder containing dicom files and/or folders|' ...
               ' Compressed file containing data|' ...
               ' Dicom files|' ...
               ' Philips PAR/REC files|' ...
               ' AFNI HEAD/BRIK files'],...
    'Position', [88 254 320 30], 'TooltipString', str, ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''srcType'');');

str = 'Enter or browse the data source according to the source type';
uicontrol('Style', 'text', 'Position', [10 210 90 30], ...
    'Tag', 'srcTxt', 'FontSize', 9, 'HorizontalAlignment', 'left', ...
    'String', 'Dicom folder', 'Background', clr, 'TooltipString', str);
uicontrol('Style', 'edit', 'Position', [88 220 296 24],'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'Background', 'white', 'Tag', 'src', ...
    'TooltipString', str, ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''set_src'');');
uicontrol('Style', 'pushbutton', 'Position', [384 221 24 22], ...
    'FontSize', 9, 'String', '...', ...
    'TooltipString', 'Browse dicom source', ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''srcDialog'');');

str = 'Enter or browse a folder to save result files';
uicontrol('Style', 'text', 'Position', [10 170 90 30], ...
    'FontSize', 9, 'HorizontalAlignment', 'left', ...
    'String', 'Result folder', 'Background', clr, 'TooltipString', str);
uicontrol('Style', 'edit', 'Position', [88 180 296 24], 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'Background', 'white', ...
    'Tag', 'dst', 'TooltipString', str, ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''set_dst'');');
uicontrol('Style', 'pushbutton', 'Position', [384 181 24 22], ...
    'FontSize', 9, 'String', '...', ...
    'TooltipString', 'Browse result folder', ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''dstDialog'');');

str = 'Choose output file format';
uicontrol('Style', 'text', 'Position', [10 132 90 30], 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'String', 'Output format', ...
    'Background', clr, 'TooltipString', str);

uicontrol('Style', 'popup', 'Background', 'white', 'Tag', 'rstFmt', ...
    'Value', 1, 'Position', [88 136 100 30], 'TooltipString', str, ...
    'String', ' .nii| .hdr/.img');

str = 'Compress into .gz files';
uicontrol('Style', 'checkbox', 'Position', [220 140 112 30], 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'String', 'Compress', ...
    'Background', clr, 'TooltipString', str, 'Tag', 'gzip');

str = 'Save one file for each volume (SPM style)';
uicontrol('Style', 'checkbox', 'Position', [330 140 72 30], 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'String', 'SPM 3D', ...
    'Background', clr, 'TooltipString', str, 'Tag', 'rst3D', ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''SPMStyle'');');
           
str = 'Choose the way to deal with SIEMENS MoCo series';
uicontrol('Style', 'text', 'Position', [10 90 90 30], 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'String', 'MoCoSeries', ...
    'Background', clr, 'TooltipString', str);
uicontrol('Style', 'popup', 'Background', 'white', 'Tag', 'mocoOpt', ...
     'Position', [88 94 320 30], 'Value', 2, 'Position', [88 94 320 30], ...
     'TooltipString', str, ...
     'String', [' Use both original and MoCo series|' ...
                ' Ignore MoCo series if both present|' ...
                ' Ignore original series if both present']);

str = ['Enter subject ID in format of Smith^John or Smith (if empty ' ...
       'first name) only if the data contains more than one subjects'];
uicontrol('Style', 'popup', 'Background', clr, 'Tag', 'subjPop', ...
    'String', ' Source is for a single subject| Convert only for subject', ...
    'Position', [10 50 210 30], ...
    'TooltipString', 'Select most likely subject information', ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''subjPop'');');
uicontrol('Style', 'edit', 'Position', [224 57 184 24], 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'Background', 'white', ...
    'Tag', 'subjName', 'TooltipString', str);

uicontrol('Style', 'pushbutton', 'Position', [104 10 200 30], ...
    'FontSize', 9, 'String', 'Start conversion', 'Tag', 'convert', ...
    'TooltipString', 'Dicom source and Result folder needed before start', ...
    'Callback', 'dicm2nii(''dicm2nii_gui_cb'',''do_convert'');');

hs = guihandles(fh); % get handles
guidata(fh, hs); % store handles
set(fh, 'HandleVisibility', 'callback'); % protect from command line

fname = [fileparts(which(mfilename)) '/dicm2nii_gui_para.mat'];
if exist(fname , 'file')
    para = load(fname); para = para.para;
    fn = fieldnames(para);
    for i = 1:length(fn)
        if strcmpi(get(hs.(fn{i}), 'Style'), 'edit')
            set(hs.(fn{i}), 'String', para.(fn{i}));
        else
            set(hs.(fn{i}), 'Value', para.(fn{i}));
        end
        dicm2nii_gui(fn{i}, hs);
    end
end
return;

%% subfunction: return phase positive in image space and fps_bits
function [phPos, fps_bits] = phaseDirection(s)
iPhase = 2; % COLUMN
foo = tryGetField(s, 'InPlanePhaseEncodingDirection', ''); % no for Philips
if strcmp(foo, 'ROW'), iPhase = 1; end
ixyz = xform_mat(s);
iPhase = ixyz(iPhase); % now can be 1/2/3 for RL/AP/IS

phPos = [];
if strncmpi(s.Manufacturer, 'SIEMENS', 7)
    phPos = csa_header(s, 'PhaseEncodingDirectionPositive'); % 0 or 1
elseif strncmpi(s.Manufacturer, 'GE', 2)
    fld = 'ProtocolDataBlock';
    if isfield(s, fld) && isfield(s.(fld), 'VIEWORDER')
        phPos = s.(fld).VIEWORDER == 1; % 1 == bottom_up
    end
elseif strncmpi(s.Manufacturer, 'Philips', 7) % no InPlanePhaseEncodingDirection
    fld = 'MRStackPreparationDirection';
    if isfield(s, 'Stack') && isfield(s.Stack.Item_1, fld)
        d = s.Stack.Item_1.(fld);
        iPhase = strfind('LRAPSIFH', d(1));
        iPhase = ceil(iPhase/2); 
        if iPhase>3, iPhase = 3; end % 1/2/3 for RL/AP/IS
        if any(d(1) == 'LPHS'), phPos = false;
        elseif any(d(1) == 'RAFI'), phPos = true;
        end
    end
end

fps_bits = [1 4 16]; % 4 for phase_dim
if iPhase == ixyz(1), fps_bits = [4 1 16]; end

%% subfunction: extract useful fields for multiframe dicom
function s = multiFrameFields(s)
if any(~isfield(s, {'SharedFunctionalGroupsSequence' ...
        'PerFrameFunctionalGroupsSequence'}))
    return; % do nothing
end
s1 = s.SharedFunctionalGroupsSequence.Item_1;
s2 = s.PerFrameFunctionalGroupsSequence.Item_1;

fld = 'EffectiveEchoTime'; n1 = 'MREchoSequence'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1; val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
end
if ~isfield(s, 'EchoTime') && ~isempty(val), s.EchoTime = val; end
if ~isfield(s, 'EchoTime') && isfield(s, 'EchoTimeDisplay')
	s.EchoTime = s.EchoTimeDisplay;
end

n1 = 'MRTimingAndRelatedParametersSequence';
fld = 'RepetitionTime'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'PixelSpacing'; n1 = 'PixelMeasuresSequence'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'SpacingBetweenSlices';  val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'SliceThickness'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'RescaleIntercept'; n1 = 'PixelValueTransformationSequence'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'RescaleSlope'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'ImageOrientationPatient'; n1 = 'PlaneOrientationSequence'; val = [];
if isfield(s1, n1)
    a = s1.(n1).Item_1;  val = tryGetField(a, fld);
elseif isfield(s2, n1)
    a = s2.(n1).Item_1;  val = tryGetField(a, fld);
end
if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end

fld = 'ImagePositionPatient'; n1 = 'PlanePositionSequence';
if isfield(s2, n1)
    a = s2.(n1); val = tryGetField(a.Item_1, fld);
    if ~isfield(s, fld) && ~isempty(val), s.(fld) = val; end
end

s2 = s.PerFrameFunctionalGroupsSequence;
n1 = fieldnames(s2);
if length(n1) < 2, return; end % in case of only 1 frame
s2 = s2.(n1{end}); % last frame

% check ImageOrientationPatient consistency for 1st and last frame only
fld = 'ImageOrientationPatient'; n1 = 'PlaneOrientationSequence';
if isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if ~isempty(val)
        try
            if any(abs(val-s.ImageOrientationPatient) > 1e-4)
                s = rmfield(s, 'ImageOrientationPatient');
                return; % inconsistent orientation, remove the field
            end
        end
    end
end

% GE data may need this to get LocationsInAcquisition
fld = 'DimensionIndexValues'; n1 = 'FrameContentSequence';
if isfield(s2, n1) && ~isfield(s, 'LocationsInAcquisition')
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if numel(val)>1, s.LocationsInAcquisition = val(2); end
end

fld = 'ImagePositionPatient'; n1 = 'PlanePositionSequence';
if isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if ~isempty(val), s.LastFile.(fld) = val; end
end

fld = 'ComplexImageComponent'; n1 = 'MRImageFrameTypeSequence';
if isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if ~isempty(val), s.LastFile.(fld) = val; end
end

fld = 'RescaleIntercept'; n1 = 'PixelValueTransformationSequence';
if isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if ~isempty(val), s.LastFile.(fld) = val; end
end
fld = 'RescaleSlope';
if isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if ~isempty(val), s.LastFile.(fld) = val; end
end

fld = 'ImagePositionPatient'; n1 = 'PlanePositionSequence';
s2 = s.PerFrameFunctionalGroupsSequence.Item_2;
if isfield(s2, n1)
    a = s2.(n1).Item_1; val = tryGetField(a, fld);
    if isfield(s, fld) && all(abs(s.(fld)-val)<1e-4)
        s.Dim3IsVolume = true;
    end
end

%% subfunction: split nii into mag and phase for Philips single file
function [nii, niiP] = split_philips_phase(nii, s)
niiP = [];
if ~strcmp(tryGetField(s, 'ComplexImageComponent', ''), 'MIXED') ... % multiframe
        && (~isfield(s, 'VolumeIsPhase') || ... 
            all(s.VolumeIsPhase) || ~any(s.VolumeIsPhase)) % not MIXED
    return;
end

if ~isfield(s, 'VolumeIsPhase') % PAR file and single-frame file have this
    dim = nii.hdr.dim(4:5);
    if tryGetField(s, 'Dim3IsVolume'), iFrames = 1:dim(2);
    else iFrames = 1:dim(1):dim(1)*dim(2);
    end
    flds = {'PerFrameFunctionalGroupsSequence' ...
        'MRImageFrameTypeSequence' 'ComplexImageComponent'};
    if dim(2) == 2 % 2 volumes, no need to re-read ComplexImageComponent
        iFrames(2) = dim(1)*dim(2); % use last frame
        s1.(flds{1}) = s.(flds{1});        
    else
        dict = dicm_dict(s.Manufacturer, flds);
        s1 = dicm_hdr(s.Filename, dict, iFrames);
    end
    s.VolumeIsPhase = false(dim(2), 1);
    for i = 1:dim(2)
        Item = sprintf('Item_%g', iFrames(i));
        foo = s1.(flds{1}).(Item).(flds{2}).Item_1.(flds{3});
        s.VolumeIsPhase(i) = strcmp(foo, 'PHASE');
    end
end

niiP = nii;
niiP.img = nii.img(:,:,:,s.VolumeIsPhase);
n = sum(s.VolumeIsPhase);
niiP.hdr.dim(5) = n; % may be 1 always
niiP.hdr.dim(1) = 3 + (n>1);

nii.img(:,:,:,s.VolumeIsPhase) = []; % now only mag
n = sum(~s.VolumeIsPhase);
nii.hdr.dim(5) = n; % may be 1 always
nii.hdr.dim(1) = 3 + (n>1);

% undo scale for 2nd set img if it was applied to img in set_nii_header
if (nii.hdr.scl_inter==0) && (nii.hdr.scl_slope==1) && ...
        (tryGetfield(s, 'RescaleIntercept') ~=0 ) && ...
        (tryGetfield(s, 'RescaleSlope') ~= 1)
    if s.VolumeIsPhase(1)
        nii.img = (nii.img - s.RescaleIntercept) / s.RescaleSlope;
        nii.hdr.scl_inter = s.LastFile.RescaleIntercept;
        nii.hdr.scl_slope = s.LastFile.RescaleSlope;
    else
        niiP.img = (niiP.img - s.RescaleIntercept) / s.RescaleSlope;
        niiP.hdr.scl_inter = s.LastFile.RescaleIntercept;
        niiP.hdr.scl_slope = s.LastFile.RescaleSlope;
    end
end

%% Write error info to a file in case user ignores Command Window output
function errorLog(errInfo)
if isempty(errInfo), return; end
fprintf(2, ' %s\n', errInfo); % red text in Command Window
global dcm2nii_errFileName;
fid = fopen(dcm2nii_errFileName, 'a');
fseek(fid, 0, -1); 
fprintf(fid, '%s\n', errInfo);
fclose(fid);

%% Clean up matlabtool related stuff: not needed for later Matlab versions
function closeMatlabpool(obj, event) %#ok 
if matlabpool('size')>0, matlabpool close; end
try stop(obj); delete(obj); end

%% Get the last date string in history
function dStr = reviseDate
persistent dateStr
if ~isempty(dateStr), dStr = dateStr; return; end
dStr = '141229?';
fid = fopen(which(mfilename));
if fid<1, return; end
str = fread(fid, '*char')';
fclose(fid);
ind = strfind(str, '% End of history. Don''t edit this line!');
if isempty(ind), return; end
ind = ind(1);
ret = str(ind-1); % new line char: \r or \n
str = str(max(1, ind-500):ind+2); % go back several lines
ind = strfind(str, [ret '% ']); % new line with % and space
for i = 1:numel(ind)-1
    ln = str(ind(i)+3 : ind(i+1)-1);
    if length(ln)>5 && all(isstrprop(ln(1:6), 'digit'))
        dStr = ln(1:6);
    end 
end
dateStr = dStr;
