function Template_matching_cortex_WU120(data_info, outDir, dmat, xdist, FSUorNU)

% 
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
%   FSUorNU: one of either 'FSU' or 'NU'; sets up institution-specific paths/parameters
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
if ~exist('FSUorNU', 'var')
    error('Must specify ''FSU'' or ''NU'' as input to set up correct paths');
end

if size(dmat) ~= [59412, 59412]
    dmat = dmat(1:59412, 1:59412);
    fprintf('Assuming distance matrix includes noncortical structures; using first 59412 grayordinates...\n')
end

dmat_thresh = double(logical(dmat > xdist));
clear dmat


%%%%%%%%%%%%%%%% add paths to Gratton Lab scripts %%%%%%%%%%%%%%%%
%addpath(genpath('/projects/b1081/Scripts'));

% set up/threshold network templates
if strcmp(FSUorNU, 'NU')
    load('/projects/b1081/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); %WU-120 consensus templates and network info
elseif strcmp(FSUorNU, 'FSU')
    load('/gpfs/research/grattonlab/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); %WU-120 consensus templates and network info
end
templates = templates(1:59412,:)';
template_values_sorted = sort(templates(:), 'descend');
threshval= template_values_sorted(round(numel(template_values_sorted) .* 0.05));
threshtemplates= templates >= threshval;
clear allTvals allTvals_sorted templates threshval

for s=1:length(subjects)
    
    fprintf('Template-matching sub-%s...\n',subjects{s})
    fprintf('   Loading subject info...\n')
    
    df = readtable(datafiles{s}); %reads into a table structure, with datafile top row as labels
    numdatas=size(df.sub,1); %number of datasets to analyses (subs X sessions)

    for i=1:numdatas

        %experiment subject IDv
        %QC(i).subjectID = df.sub{i}; 
        if isa(df.sub(i),'cell')
            QC(i).subjectID = df.sub{i}; % the more expected case
        elseif isa(df.sub(i),'double')  %to account for subject numbers that are all numeric
            QC(i).subjectID = num2str(df.sub(i)); %change to string to work with rest of code
        else
            error('can not recognize subject data type')
        end
        
        if ~strcmp(subjects{s}, QC(i).subjectID)
            error('Make sure subjects column aligns with all rows in Datalist.xslx file')
        end
        QC(i).session = df.sess(i); %session ID
        QC(i).condition = df.task{i}; %condition type (rest or name of task)
        QC(i).TR = df.TR(i,1); %TR (in seconds)
        QC(i).TRskip = df.dropFr(i); %number of frames to skip
        QC(i).topDir = df.topDir{i}; %initial part of directory structure
        QC(i).dataFolder = df.dataFolder{i}; % folder for data inputs (assume BIDS style organization otherwise)
        QC(i).confoundsFolder = df.confoundsFolder{i}; % folder for confound inputs (assume BIDS organization)
        QC(i).FDtype = df.FDtype{i,1}; %use FD or fFD for tmask, etc?
        QC(i).runs = str2double(regexp(df.runs{i},',','split'))'; % get runs, converting to numerical array (other orientiation since that's what's expected)    
        
        % to address potential residuals feild (to indicate residuals for task FC)
        % if ismember('residuals',df.Properties.VariableNames)
        %     QC(       i).residuals = df.residuals(i);
        % else
        %     QC(i).residuals = 0; % assume these are not residuals if the field doesn't exist
        % end
    end
    
    catData = [];
    
    fprintf('   Concatenating BOLD (CIFTI) rest data...\n')
    
    for i = 1:numdatas
	
	%for now... skip non-rest data	
	if ~strcmp(QC(i).condition,'rest')
	    continue
	end
	
        for r = 1:length(QC(i).runs) % loop through runs (may not be continuous)
            
            conf_fstring1 = sprintf('%s/%s/fmriprep/sub-%s/ses-%d/func/',QC(i).topDir,QC(i).confoundsFolder,QC(i).subjectID,QC(i).session);
            all_fstring2 = sprintf('sub-%s_ses-%d_task-%s_run-%d',QC(i).subjectID,QC(i).session,QC(i).condition,QC(i).runs(r));
            cifti_fstring = sprintf('%s/sub-%s/ses-%d/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',cifti_dirs{s},QC(i).subjectID,QC(i).session,QC(i).subjectID,QC(i).session,QC(i).runs(r));
            
            %
            %   ***residuals -- not yet set up***
            %
            % if QC(i).residuals == 0 % the typical case
            %     bolddata_fname = [data_fstring1 all_fstring2 '_space-' space '_' res '_desc-preproc_bold.nii.gz'];
            % else %if these are task FC that have been residualized
            %     bolddata_fname = [data_fstring1 all_fstring2 '_space-' space '_' res '_desc-preproc_bold_residuals.nii.gz'];
            % end
            
            tmask_fname = [conf_fstring1 'FD_outputs/' all_fstring2 '_desc-tmask_' QC(i).FDtype '.txt']; %assume this is in confounds folder
            boldmot_folder{i,r} = [conf_fstring1 'FD_outputs']; % in this case, just give path/start so I can load different versions

            if ~exist(boldmot_folder{i,r})
                error(['FD folder does not exist for: ' boldmot_folder{i,r}]);
            end
            
            run_tmask = table2array(readtable(tmask_fname));
            run_data = ft_read_cifti_mod(cifti_fstring);

            catData = [catData run_data.data(1:59412, run_tmask>0)];
            
            if i==1 && r==1
                template_cifti = run_data; template_cifti.data=[];
            end
            
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
    dice_match_fname = sprintf('%s/sub-%s_dice_to_templates.mat',outDir,QC(i).subjectID);
    save(dice_match_fname, 'dice_to_templates','-v7.3');
    dice_map_fname = sprintf('%s/sub-%s_dice_WTA_map_kden0.05.dtseries.nii',outDir,QC(i).subjectID);
    template_cifti.data = dice_subject_index;
        template_cifti.time=1; template_cifti.hdr.dim(6)=1; template_cifti.hdr.dim(7)=59412; template_cifti.brainstructure=template_cifti.brainstructure(1:64984,:);
        template_cifti.brainstructurelabel={'CORTEX_LEFT','CORTEX_RIGHT'}; template_cifti.pos=template_cifti.pos(1:64984,:);
    ft_write_cifti_mod(dice_map_fname, template_cifti)
    clear template_cifti
end

end