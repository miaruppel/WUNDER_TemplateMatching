function smooth_volume_222(datafile,FCdir,blurkernel)
%
% this function applies a Gaussian blur to all FCprocessed BOLD runs
% specified in the input Datalist, and output a volume file into the
% specified FCprocess directory (same FCprocess directory with unsmoothed data)
%
% "blurkernel" input is optional; IF NOT SPECIFIED AS THIRD INPUT, defaults
% to 4 (millimeters)
% (4 is standard for data in 222 resolution. For data in 333 resolution, use 6)
%
% EXAMPLEs: 
% smooth_volume_222('EXAMPLESUB_DATALIST.xlsx','/projects/b1081/iNetworks/Nifti/derivatives/preproc_FCProc/')
% smooth_volume_222('EXAMPLESUB_DATALIST.xlsx','/projects/b1081/iNetworks/Nifti/derivatives/preproc_FCProc/',6)
%
% 12.10.2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SPECIFY BLUR KERNEL (in mm); default is 4 (standard for data in 222)
if nargin == 2
    blurkernel = 4;
end
fprintf(['\nSmoothing data with blur kernel of ' num2str(blurkernel) 'mm\n']);


%% load in sub/session info
df = readtable(datafile); %reads into a table structure, with datafile top row as labels
numdatas=size(df.sub,1); %number of datasets to analyses (subs X sessions)

% organize for easier use
for i = 1:numdatas
    
    % reformat sub and run info as needed
    if isa(df.sub(i),'cell')
        subInfo(i).subjectID = df.sub{i}; % the more expected case
    elseif isa(df.sub(i),'double')  %to account for subject numbers that are all numeric
        subInfo(i).subjectID = num2str(df.sub(i)); %change to string to work with rest of code
    else
        error('can not recognize subject data type')
    end
    
    subInfo(i).session = df.sess(i); %session ID
    subInfo(i).condition = df.task{i}; %condition type (rest or name of task)
    subInfo(i).TR = df.TR(i,1); %TR (in seconds)
    subInfo(i).TRskip = df.dropFr(i); %number of frames to skip
    subInfo(i).topDir = df.topDir{i}; %initial part of directory structure
    subInfo(i).dataFolder = df.dataFolder{i}; % folder for data inputs (assume BIDS style organization otherwise)
    subInfo(i).confoundsFolder = df.confoundsFolder{i}; % folder for confound inputs (assume BIDS organization)
    subInfo(i).FDtype = df.FDtype{i,1}; %use FD or fFD for tmask, etc?
    subInfo(i).runs = str2double(regexp(df.runs{i},',','split'))'; % get runs, converting to numerical array (other orientiation since that's what's expected
    %subInfo(i).space = space; %assuming same space for now... insert check
    %subInfo(i).res = res; %assuming same res for now.. insert check
end

%% Loop through data, extract timecourses
for i = 1:numdatas
    
    fprintf('Subject %s, session %d\n',subInfo(i).subjectID,subInfo(i).session);
    
    
    for j = 1:length(subInfo(i).runs)
        
        fprintf('Run: %d\n',subInfo(i).runs(j));
        
        % set FCprocessed and output (smoothed) filenames:
        FCprocFile = sprintf('%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-%s_run-%d_fmriprep_zmdt_resid_ntrpl_bpss_zmdt.nii.gz',...
            FCdir,subInfo(i).subjectID,subInfo(i).session,...
            subInfo(i).subjectID,subInfo(i).session,subInfo(i).condition,subInfo(i).runs(j));
        smoothProcFile = sprintf(['%s/sub-%s/ses-%d/func/sub-%s_ses-%d_task-%s_run-%d_fmriprep_zmdt_resid_ntrpl_bpss_zmdt_g' num2str(blurkernel) 'mm.nii.gz'],...
            FCdir,subInfo(i).subjectID,subInfo(i).session,...
            subInfo(i).subjectID,subInfo(i).session,subInfo(i).condition,subInfo(i).runs(j));
        
        % run AFNI command to smooth data
        system(['module load singularity; singularity exec -B /projects/b1081:/projects/b1081 /projects/b1081/singularity_images/afni_latest.sif 3dmerge -1blur_fwhm ' num2str(blurkernel) ' -doall -prefix ' smoothProcFile ' ' FCprocFile ' -overwrite']);	

    end
end

end

