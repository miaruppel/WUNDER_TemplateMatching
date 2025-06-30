function make_fs_masks(subject, fmriprepTopDir, make_fs_make_params_file)

%%%%%%%%%%%%%%%%%%%%%%
% This makes several erosions of WM masks, including no erosion.
% JDP 8/14/12, modified from TOL script; modified for NU/Matlab Oct. 2020
%
% The outputs (e.g., 'sub-INET003_space-MNI152NLin6Asym_label-WM_probseg_0.9mask_res-2_ero3.nii.gz')
% are written into the subject's overall 'anat' folder.
%
% First input is subject ID (e.g., 'INET003'; assumes BIDS structure); second is fmriprep
% directory (e.g., '/projects/b1081/iNetworks/Nifti/derivatives/preproc_fmriprep-20.2.0');
% third is probabilistic threshold to start erosions
% also assumes fmriprep directory structure as our standard
%

%%%%%% change parameters if desired %%%%%%
space = 'MNI152NLin6Asym';
voxdim = '2'; %voxel size
eroiterwm = 4; %number of erosions to perform
WMprobseg_thresh = 0.9;
%------------
WMprobseg = ['sub-' subject '_space-' space '_label-WM_probseg.nii.gz'];
WMmaskname = ['sub-' subject '_space-' space '_label-WM_probseg_' num2str(WMprobseg_thresh) 'mask.nii.gz'];
anat_dir= [fmriprepTopDir '/fmriprep/sub-' subject '/anat/'];


% Find params file and read in system-specific variables
[paramspath,paramsname,paramsextension] = fileparts(make_fs_make_params_file);
origpath = pwd;
if ~isempty(paramspath)
    cd(paramspath)
end
%Load parameters
params = feval(paramsname);
varnames = fieldnames(params);
for i = 1:length(varnames)
    evalc([varnames{i} ' = params.' varnames{i}]);
end
clear varnames params
cd(origpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(anat_dir, 'dir')
    cd(anat_dir);
else
    disp('No overall anat directory; assuming subject has 1 scan session with a T1')
    
    %find session with anat dir
    sessions = dir([fmriprepTopDir '/fmriprep/sub-' subject '/ses-*']); sessions = sessions([sessions(:).isdir]);
    sessions = {sessions.name};
    for s=1:size(sessions,2)
        if exist([fmriprepTopDir '/fmriprep/sub-' subject '/' sessions{s} '/anat'],'dir')
            anat_ses = sessions{s};
            break
        end
    end
    if ~exist('anat_ses', 'var')
        error('No session found with an anat directory')
    end
    
    %try to locate existing T1; want the one in native space
    T1w_filenames = dir([fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/*desc-preproc_T1w.nii.gz']); T1w_filenames = {T1w_filenames.name};
    T1w_filenames = T1w_filenames(cellfun(@(s)isempty(regexp(s,'space')),T1w_filenames)); %eliminates any T1s in atlas space
    
    if size(T1w_filenames,2) == 0
        error('No T1 images found in session anat folder; check T1s and retry')
    elseif size(T1w_filenames,2) > 1
        error('More than 1 usable T1 image found in session anat folder; consider linking files manually')
    else
    % make anat dir, link necessary files to those in session with anat folder; 
    % masks (probseg files) should be in res-02 by default from fmriprep
        disp(['T1 in native space identified as: ' T1w_filenames{1}])
        disp('    Continuing...')
        mkdir(anat_dir); cd(anat_dir);
        system(['ln -s ' fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/' T1w_filenames{1} ' sub-' subject '_desc-preproc_T1w.nii.gz']);
        h5_filename = dir([fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/*from-T1w_to-' space '_mode-image_xfm.h5']);
        h5_filename = h5_filename(1).name;
        system(['ln -s ' fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/' h5_filename ' sub-' subject '_from-T1w_to-' space '_mode-image_xfm.h5']);
        
        %this should point to the WM_probseg in native space (i.e. not MNI152NLin6Asym)
        WM_filename = dir([fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/*WM_probseg.nii.gz']);
        WM_filename = WM_filename(1).name;
        system(['ln -s ' fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/' WM_filename ' sub-' subject '_label-WM_probseg.nii.gz']);
        
        %other masks only need to be linked in their new space/resolution for use in FCprocess
        GM_filename = dir([fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/*_space-' space '_res-' voxdim '_*GM_probseg.nii.gz']);
        GM_filename = GM_filename(1).name;
        system(['ln -s ' fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/' GM_filename ' sub-' subject '_space-' space '_res-' voxdim '_label-GM_probseg.nii.gz']);
        CSF_filename = dir([fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/*_space-' space '_res-' voxdim '_*CSF_probseg.nii.gz']);
        CSF_filename = CSF_filename(1).name;
        system(['ln -s ' fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/' CSF_filename ' sub-' subject '_space-' space '_res-' voxdim '_label-CSF_probseg.nii.gz']);
        brainmask_filename = dir([fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/*_space-' space '_res-' voxdim '_*brain_mask.nii.gz']);
        brainmask_filename = brainmask_filename(1).name;
        system(['ln -s ' fmriprepTopDir '/fmriprep/sub-' subject '/' anat_ses '/anat/' brainmask_filename ' sub-' subject '_space-' space '_res-' voxdim '_desc-brain_mask.nii.gz']);
    end
end


%%% first, use ANTs to warp T1 and WMprobseg images into MNI in 111 space, if they don't yet exist %%%
T1_templateLoc = [templateflow_dir '/tpl-' space '/tpl-' space '_res-01_T1w.nii.gz'];
inNames = {['sub-' subject '_desc-preproc_T1w.nii.gz'], ['sub-' subject '_label-WM_probseg.nii.gz']};
outNames = {['sub-' subject '_space-' space '_desc-preproc_T1w.nii.gz'], WMprobseg};

for tform = 1:length(inNames)
    if ~exist(outNames{tform}, 'file')
        system([singularity_cmd_start 'singularity exec -B ' project_dir ':' project_dir ',' anat_dir ' ' ants_sif ' antsApplyTransforms --verbose -i ' inNames{tform} ' -o ' outNames{tform} ' -r ' T1_templateLoc ' -t sub-' subject '_from-T1w_to-' space '_mode-image_xfm.h5']);
    end
end


%%% threshold at WMprobseg_thresh and binarize %%%
system([fsl_cmd_start 'fslmaths ' WMprobseg  ' -thr ' num2str(WMprobseg_thresh) ' -bin ' WMmaskname]);


%%% erode cerebral WM mask to avoid possible gray matter contamination %%%
iter = 0;
%system(['module load singularity; singularity run -B ' anat_dir ' /projects/b1081/singularity_images/afni_latest.sif 3dresample -dxyz ' voxdim ' ' voxdim ' ' voxdim ' -prefix ' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz -input ' WMmaskname]);
system([singularity_cmd_start 'singularity exec -B ' project_dir ':' project_dir ',' anat_dir ' ' afni_sif ' 3dresample -dxyz ' voxdim ' ' voxdim ' ' voxdim ' -prefix ' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz -input ' WMmaskname]);
system([fsl_cmd_start 'fslmaths ' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz -bin ' WMmaskname(1:end-7) '_res-' voxdim '_ero0.nii.gz']);

iter = 1;
while iter <= eroiterwm
    system([fsl_cmd_start 'fslmaths ' WMmaskname ' -kernel 3D -ero ' WMmaskname]);
    	
    %system(['module load singularity; singularity exec -B /projects/b1081:/projects/b1081,' anat_dir ' /projects/b1081/singularity_images/afni_latest.sif 3dresample -dxyz ' voxdim ' ' voxdim ' ' voxdim ' -prefix ' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz -input ' WMmaskname ' -overwrite']);
    system([singularity_cmd_start 'singularity exec -B ' project_dir ':' project_dir ',' anat_dir ' ' afni_sif ' 3dresample -dxyz ' voxdim ' ' voxdim ' ' voxdim ' -prefix ' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz -input ' WMmaskname ' -overwrite']);

    system([fsl_cmd_start 'fslmaths ' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz -bin ' WMmaskname(1:end-7) '_res-' voxdim '_ero' num2str(iter) '.nii.gz']);
    iter = iter + 1;
end

%%% remove unnecessary files %%%
system(['rm ' anat_dir '/' WMmaskname ' ' anat_dir '/' WMmaskname(1:end-7) '_res-' voxdim '.nii.gz']);
cd(origpath)

end
