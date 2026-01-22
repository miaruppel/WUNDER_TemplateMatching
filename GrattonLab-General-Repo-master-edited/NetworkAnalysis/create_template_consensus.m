% create new template
% HCP dtseries files to template_consensus.mat 

network_names = {'DMN', 'Vis', 'FP', 'DAN', 'Sal', 'CO', 'SMd', 'SMl', 'Aud', 'PMN', 'PON'}; % network names
IDs = {1, 2, 3, 5, 8, 9, 10, 11, 12, 15, 16}; % unique IDs for each network
templates = zeros(59412, length(network_names), 'single'); % initialize templates matrix

% read in dtseries files 
dtseries_folder =  '/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/TemplateMatching/HCP_CIFTI_templates';

for i = 1:length(network_names)
    % filepath for the current dtseries file
    dtseries_file = fullfile(dtseries_folder, sprintf('HCP384_%s_template.dtseries.nii', network_names{i}));
    
    % load the dtseries file
    cifti_data = ft_read_cifti_mod(dtseries_file); % Replace with appropriate reader if not using FieldTrip

    % extract data and populate templates matrix
    % assuming data dimension is [vertices x 1]
    templates(:, i) = single(cifti_data.data); % Ensure data is in single precision
end
