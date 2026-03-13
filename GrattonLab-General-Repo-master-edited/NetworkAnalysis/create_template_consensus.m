% create new template
% ABCD dscalar files to template_consensus.mat 

network_names = {'DMN', 'Vis', 'FP', 'DAN', 'VAN', 'Sal', 'CO', 'SMd', 'SMl', 'Aud', 'Tpole', 'MTL', 'PMN', 'PON'}; % network names
IDs = {1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}; % unique IDs for each network

% dscalar files path
dscalar_folder =  '/data/wheelock/data1/parcellations/Hermosillo_2024_NatureNeuroscience/probability_maps';

% Read first file to get vertex count
tmp = ft_read_cifti_mod(fullfile(dscalar_folder, ...
    'ABCD_GRP1_singlenet_probability_n2988_Aud_network.dscalar.nii'));
nVert = size(tmp.data,1);
templates = zeros(nVert, length(network_names), 'single');

for i = 1:length(network_names)
    % filepath for the current dscalar file
    fname = sprintf('ABCD_GRP1_singlenet_probability_n2988_%s_network.dscalar.nii', ...
                    network_names{i});
    fpath = fullfile(dscalar_folder, fname);
    
    % load the dtseries file
    cifti_data = ft_read_cifti_mod(fpath); % Replace with appropriate reader if not using FieldTrip

    % extract data and populate templates matrix
    % assuming data dimension is [vertices x 1]
    templates(:, i) = single(cifti_data.data); % Ensure data is in single precision
end
