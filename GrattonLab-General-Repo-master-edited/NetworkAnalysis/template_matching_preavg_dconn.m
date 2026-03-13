% template matching with pre-averaged dconns!

% add paths to toolboxes 
addpath('/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/NetworkAnalysis/cifti-matlab-master');
addpath('/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/NetworkAnalysis/functions');
addpath('/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/NetworkAnalysis/fieldtrip-20240731');

% load in template (WU-120 or HCP or ABCD)
%load('/data/cn/data1/scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'); %WU-120 consensus templates and network info
%load('/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/TemplateMatching/HCP_CIFTI_templates/Templates_consensus.mat'); % generated HCP templates
load('/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/ABCD_template/Hermosillo_ABCD_template_n141_zscored.mat'); % Robert's ABCD template (includes SCAN as 15th network)
templates = templates(1:59412,:)';
template_values_sorted = sort(templates(:), 'descend');
threshval= template_values_sorted(round(numel(template_values_sorted) .* 0.05));
threshtemplates= templates >= threshval;
clear allTvals allTvals_sorted templates threshval

% load in averaged dconn data
cifti_path = '/data/shimony/shimony2/wunder/wunder_caf_II/ABCD_DCAN_docker/docker_output_BP_filter_and_cleaning_2020_05_17/dconns/N52_termcontrols_group_avg.dconn.nii';
avg_dconn = ft_read_cifti_mod(cifti_path);
dtseries_path = '/data/smyser/smyser4/wunder/wunder_caf_III/recon_docker/docker_output/sub-caf067f/ses-None/files/MNINonLinear/Results/task-rest_DCANBOLDProc_v4.0.0_Atlas_smooth2.55.dtseries.nii';
template_cifti = ft_read_cifti_mod(dtseries_path);
%template_cifti = avg_dconn; 
template_cifti.data = [];

data = avg_dconn.data(1:59412, 1:59412);
data = single(data);
data_vals = data(triu(true(59412),1));
sorted_vals = sort(data_vals, 'descend');
subject_thresh = sorted_vals(round(0.05 * numel(sorted_vals)));
corr_mat_thresh = data >= subject_thresh;
clear data data_vals sorted_vals

fprintf('   Calculating dice overlap to templates...\n')
%%% compute dice overlap b/w each vertex and template %%%
dice_to_templates= zeros(length(IDs),59412);
for t=1:length(IDs)        
    big_cifti_template_mat_full= repmat(threshtemplates(t,:),59412,1);
    dice_to_templates(t,:)= (sum((big_cifti_template_mat_full & corr_mat_thresh),2)*2) ./ (sum(big_cifti_template_mat_full,2) + sum(corr_mat_thresh,2));
end
clear corr_mat_thresh

%%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
[x, dice_subject_index] = max(dice_to_templates',[],2);
dice_subject_index = (IDs(dice_subject_index))';
dice_subject_index(x==0)=0;

%%% write out results %%%
outDir = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/output_910_data/termcontrols/ABCD_template';
QC.subjectID = 'termcontrols_avgdconn_ABCD';

dice_match_fname = sprintf('%s/sub-%s_0.2FD_dice_to_templates.mat',outDir,QC.subjectID);
save(dice_match_fname, 'dice_to_templates','-v7.3');
dice_map_fname = sprintf('%s/sub-%s_0.2FD_dice_WTA_map_kden0.05.dtseries.nii',outDir,QC.subjectID);
template_cifti.data = dice_subject_index;
    template_cifti.time=1; template_cifti.hdr.dim(6)=1; template_cifti.hdr.dim(7)=59412; template_cifti.brainstructure=template_cifti.brainstructure(1:64984,:);
    template_cifti.brainstructurelabel={'CORTEX_LEFT','CORTEX_RIGHT'}; template_cifti.pos=template_cifti.pos(1:64984,:);
ft_write_cifti_mod(dice_map_fname, template_cifti)
