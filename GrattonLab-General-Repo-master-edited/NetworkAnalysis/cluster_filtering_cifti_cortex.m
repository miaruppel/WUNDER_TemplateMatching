function cluster_filtering_cifti_cortex(networks_file, clustersize_mm, L_surface, R_surface, out_dir)

%
% inputs: networks_file: path to CIFTI network file (this is set up to take a *dtseries*)
%       clustersize_mm: surface area (in millimeters); all clusters smaller than this size will become 0 in final map
%       out_dir: place to (temporarily) store intermediate files and final clustered network map
%       L_surface & R_surface: paths to surfaces on which to calculate surface area (e.g., Conte midthickness surfaces or subject-specific)
% outputs: writes out a CIFTI (dtseries) into output folder with clusters of size <[clustersize_mm] zero'd out
%
% example:
% cluster_filtering_cifti_cortex('/path/to/subject_dice_to_templates_map.dtseries.nii', 50, ...
%   '/path/to/CConte69.L.midthickness.32k_fs_LR.surf.gii', '/path/to/Conte69.L.midthickness.32k_fs_LR.surf.gii', '/path/to/output_dir')

%check/format inputs
if ~strcmp(networks_file(end-12:end), '.dtseries.nii')
    error('input networks file must be a CIFTI dtseries file')
end
[~, filename, ~] = fileparts(networks_file);

networks = ft_read_cifti_mod(networks_file); 
networks_data = networks.data;
network_nums = unique(networks.data); network_nums(network_nums==0) = [];

disp(['Clusterizing ' networks_file ' at ' num2str(clustersize_mm) 'mm...'])

%%% separate networks by cifti %%%
bin_networks = zeros(size(networks_data,1),length(network_nums));
for i=1:length(network_nums)
    bin_networks(:,i) = double(networks_data == network_nums(i));
end
tempcifti = networks;
tempcifti.data = bin_networks;
temp_filename = [out_dir '/temp_separated_by_networks.dtseries.nii'];
ft_write_cifti_mod(temp_filename, tempcifti);

%%% clusterize and size exclude %%%
value_thresh = '0'; %all values wil be ones
tempclusters_filename = [out_dir '/temp_clusters_by_network.dtseries.nii'];
find_clusters = ['wb_command -cifti-find-clusters ' temp_filename ' ' value_thresh ' ' num2str(clustersize_mm) ' 0 0 COLUMN ' tempclusters_filename ' -left-surface ' L_surface ' -right-surface ' R_surface];
system(find_clusters);             
clear filepath outputfilepath

clusters = ft_read_cifti_mod(tempclusters_filename);
clustered = networks; clustered.data = zeros(size(clustered.data));
for i = 1:length(network_nums)
    clustered.data(:,i) = double(logical(clusters.data(:,i) > 0)) * network_nums(i);
end
clustered.data = sum(clustered.data, 2);

clustered_filename = [out_dir '/' filename(1:end-9) '_clusterMinSize' num2str(clustersize_mm) '.dtseries.nii'];
ft_write_cifti_mod(clustered_filename, clustered)
system(['rm ' out_dir '/temp_separated_by_networks.dtseries.nii ' out_dir '/temp_clusters_by_network.dtseries.nii']);


end
