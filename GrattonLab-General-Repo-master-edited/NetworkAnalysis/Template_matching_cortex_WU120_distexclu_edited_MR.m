function Template_matching_cortex_WU120_distexclu_edited_MR(data_info, outDir, xdist, FSUorHCPorWU)

% Template_matching_cortex_WU120_edited_MR('data_info_preemies_lowmotion_noinjury_batch.txt', '/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/TemplateMatching/output/preemies_lowmotion_noinjury_0.2FD_batch_HCP_template', 'HCP')
% Matches each vertex to a functional network by comparing the dice coeff.
% of its thresholded seedmap to all network template seedmaps (from the
% WashU 120), and assigns the vertex to the network with the highest match.
% (assumes both template and target surfaces are the 59,412-vertex cortical surface)
%
% Inputs:
%   data_info: path to a .txt file with 3 columns: (1) a list of subject 
%   IDs (e.g., 'INET002'); (2) the path to the corresponding subject's .xlsx
%   datalist (see GrattonLab github for example datalist); (3) the path to
%   the top directory which contains subjects and their CIFTI (assumed file
%   structure of, e.g., /projects/b1081/Lifespan/postFCproc_CIFTI/sub-LS03/
%   ses-1/cifti_timeseries_normalwall/sub-LS03_XXXXXXX_smooth2.55.dtseries.nii)
%   example row of data_info.txt:   LS03    /home/amd3565/Lifespan_Datalists/lifespan_datalist_LS03.xlsx    /projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI
%
%   outDir: output destination for dtseries.nii map and Dice values .mat
%
%   FSUorNUorWU: one of either 'FSU' or 'NU'; sets up institution-specific paths/parameters
%                   WU component added 2024-08-28 (MR)
%
%   *** note: this script currently assumes runs labeled as 'rest' in the .xlsx
%       datalist-- be sure to use a datalist only including rest runs ***
%	(not currently set up for residuals either. talk to Ally if you need this)
%
% Outputs:
%   (1) a CIFTI sub-{subject}_dice_WTA_map_kden0.05.dtseries.nii, assigning
%       each vertex to one of 14 networks as coded by Power colors (see info
%       from loaded Templates_consensus.mat file), based on highest dice coeff.
%       
%   (2) a sub-{subject}_dice_to_templates.mat, which saves the dice coeff. of
%       each vertex's thresholded seedmap to each template's seedmap.
%
%
%

[subjects, datafiles, cifti_dirs] = textread(data_info, '%s%s%s');
if ~exist('FSUorHCPorWU', 'var')
    error('Must specify ''FSU'' or ''HCP'' or ''WU''as input to set up correct paths');
end

% Conte reference distance matrix
%load('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_large_uint8.mat', 'distances');
%if size(distances) ~= [59412, 59412]
%    distances = distances(1:59412, 1:59412);
%    fprintf('Assuming distance matrix includes noncortical structures; using first 59412 grayordinates...\n')
%end

%%%%%%%%%%%%%%%% add paths to Gratton Lab scripts %%%%%%%%%%%%%%%%
% for some reason, still having to specicify these in the command window for this to work 
addpath('cifti-matlab-master');
addpath('functions');
addpath('fieldtrip-20240731');


% set up/threshold network templates
if strcmp(FSUorHCPorWU, 'HCP')
    load('/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/TemplateMatching/HCP_CIFTI_templates/Templates_consensus.mat'); % generated HCP templates
elseif strcmp(FSUorHCPorWU, 'FSU')
    load('/gpfs/research/grattonlab/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); 
elseif strcmp(FSUorHCPorWU, 'WU')
    load('/data/cn/data1/scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); %WU-120 consensus templates and network info
end
templates = templates(1:59412,:)';
template_values_sorted = sort(templates(:), 'descend');
threshval= template_values_sorted(round(numel(template_values_sorted) .* 0.05));
threshtemplates= templates >= threshval;
clear allTvals allTvals_sorted templates threshval

for s=1:length(subjects)
    fprintf('Subject #%d \n', s)
    fprintf('Template-matching sub-%s...\n',subjects{s})
    fprintf('   Loading subject info...\n')
    
    df = readtable(datafiles{s},'PreserveVariableNames', true); %reads into a table structure, with datafile top row as labels
   
    QC.subjectID = df.sub{s};  % Assuming the first row should match the subject
    QC.session = df.sess(s);  % Session ID
    QC.condition = df.task{s};  % Condition type
    QC.TR = df.TR(s);  % TR in seconds
    QC.TRskip = df.dropFr(s);  % Number of frames to skip
    QC.topDir = df.topDir{s};  % Initial part of directory structure
    QC.dataFolder = df.dataFolder{s};  % Folder for data inputs
    QC.confoundsFolder = df.confoundsFolder{s};  % Folder for confound inputs
    QC.FDtype = df.FDtype{s};  % FD type for tmask, etc.
    QC.runs = str2double(regexp(df.runs{s}, ',', 'split'))';
    %QC.runs = 1;
    
    catData = [];
    
    fprintf('   Concatenating BOLD (CIFTI) rest data...\n')
    %save('debug_workspace1.mat');
    
    % load in subject-specific distance matrix 
    dist_matrix_path = sprintf('/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master/NetworkAnalysis/%s_distance_matrix.mat', QC.subjectID);
    load(dist_matrix_path, 'M');
    
    dmat_thresh = double(logical(M > xdist));
    clear M
    
  
        for r = 1:length(QC.runs) % loop through runs (may not be continuous)
            fprintf('Run #%d \n', r)
            conf_fstring1 = sprintf('%s/%s/sub-%s/ses-None/files/',QC.topDir,QC.confoundsFolder,QC.subjectID);
            all_fstring2 = sprintf('sub-%s_task-%s_run-%d',QC.subjectID,QC.condition,QC.runs(r));
            cifti_fstring = sprintf('%s/docker_output/sub-%s/ses-None/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Atlas_smooth2.55.dtseries.nii',cifti_dirs{s},QC.subjectID);
            %save('debug_workspace.mat')
            fprintf('%s \n', cifti_fstring)
            
            %   ***residuals -- not yet set up***
            %
            % if QC(i).residuals == 0 % the typical case
            %     bolddata_fname = [data_fstring1 all_fstring2 '_space-' space '_' res '_desc-preproc_bold.nii.gz'];
            % else %if these are task FC that have been residualized
            %     bolddata_fname = [data_fstring1 all_fstring2 '_space-' space '_' res '_desc-preproc_bold_residuals.nii.gz'];
            % end
            
            %mins = 20;
            %all_fstring3 = sprintf('sub-%s_tmask_%dmins',QC.subjectID,mins);
            %tmask_fname = [conf_fstring1 'DCANBOLDProc_v4.0.0/analyses_v2/motion/' all_fstring3 '.txt']; %assume this is in confounds folder

            tmask_fname = [conf_fstring1 'DCANBOLDProc_v4.0.0/analyses_v2/motion/' all_fstring2 '-tmask_0.2' QC.FDtype '.txt']; %assume this is in confounds folder
            boldmot_folder{s,r} = [conf_fstring1 'DCANBOLDProc_v4.0.0']; % in this case, just give path/start so I can load different versions

            if ~exist(boldmot_folder{s,r})
                error(['FD folder does not exist for: ' boldmot_folder{s,r}]);
            end
            
            run_tmask = table2array(readtable(tmask_fname));
            %save('debug_workspace2.mat');
            run_data = ft_read_cifti_mod(cifti_fstring);

            catData = [catData run_data.data(1:59412, run_tmask>0)];
            
            if s==1 && r==1
                template_cifti = run_data; template_cifti.data=[];
            end
            
        end
    
    %%% compute subject's correlation matrix %%%
    fprintf('   Computing vertexwise correlation matrix on data of size %d by %d...\n',size(catData,1),size(catData,2))
    corr_mat_full = paircorr_mod(catData');
    corr_mat_full = single(corr_mat_full);
    corr_vals = corr_mat_full(triu(true(59412),1));
    sorted_vals = sort(corr_vals, 'descend');
    subject_thresh = sorted_vals(round(0.05 * numel(sorted_vals)));
    corr_mat_thresh = corr_mat_full >= subject_thresh;
    corr_mat_thresh_xdist = (corr_mat_thresh & dmat_thresh);
    
    clear catData corr_mat_full corr_vals sorted_vals
    
    fprintf('   Calculating dice overlap to templates...\n')
    
    %%% compute dice overlap b/w each vertex and template %%%
    dice_to_templates= zeros(length(IDs),59412);
    for t=1:length(IDs)        
        big_cifti_template_mat_full= repmat(threshtemplates(t,:),59412,1);
        dice_to_templates(t,:)= (sum((big_cifti_template_mat_full & corr_mat_thresh_xdist),2)*2) ./ (sum(big_cifti_template_mat_full,2) + sum(corr_mat_thresh_xdist,2));
    end
    clear corr_mat_thresh corr_mat_thresh_xdist

    %%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
    [x, dice_subject_index] = max(dice_to_templates',[],2);
    dice_subject_index = (IDs(dice_subject_index))';
    dice_subject_index(x==0)=0;
    
    %%% write out results %%%
    dice_match_fname = sprintf('%s/sub-%s_0.2FD_dice_to_templates.mat',outDir,QC.subjectID);
    save(dice_match_fname, 'dice_to_templates','-v7.3');
    dice_map_fname = sprintf('%s/sub-%s_0.2FD_dice_WTA_map_kden0.05.dtseries.nii',outDir,QC.subjectID);
    template_cifti.data = dice_subject_index;
        template_cifti.time=1; template_cifti.hdr.dim(6)=1; template_cifti.hdr.dim(7)=59412; template_cifti.brainstructure=template_cifti.brainstructure(1:64984,:);
        template_cifti.brainstructurelabel={'CORTEX_LEFT','CORTEX_RIGHT'}; template_cifti.pos=template_cifti.pos(1:64984,:);
    ft_write_cifti_mod(dice_map_fname, template_cifti)
    %clear template_cifti
    
    fprintf('Finished: sub-%s \n', subjects{s})
  
end

clear template_cifti

end