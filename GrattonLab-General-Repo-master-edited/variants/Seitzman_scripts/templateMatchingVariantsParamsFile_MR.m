function params = templateMatchingVariantsParamsFile_MR
%Parameters to set for template_matching_variants.m
%
% All required files are detailed below, with example paths provided.

% List of subject timecourses in CIFTI format. The file should contain two
% columns separated by a space: subjNumber timecourseLocation
surfdatafile = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/variants/Seitzman_scripts/N14_preemies_datalist.txt'; 

% List of temporal masks for the timecourses (frame censoring). The file
% should contain two columns separated by a space: subjNumber maskLocation
tmaskfile = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/variants/Seitzman_scripts/N14_preemies_tmasklist.txt';

% List of network variants. The file should contain one column: variantLocation
variantsfile = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/variants/Seitzman_scripts/N14_preemies_variantslist.txt';

% Location of template networks. Should contain a 'templates' variable that
% is (#vertices) X (#networks) and contains unthresholded connectivity maps
% per network. Should also contain 'IDs', a 1 X (#systems) vector of 
% numeric IDs for each system.
templatesfile = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/GrattonLab-General-Repo-master-edited/variants/Seitzman_scripts/networkTemplates.mat'; 


%--------------------------------------------------------------------------

% Load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end
