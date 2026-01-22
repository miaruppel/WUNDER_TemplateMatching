% Wrapper script to loop over subjects for createSpatialCorrMap_MR.m

% Path to text file with subject dconn locations
% 9/10 data
% terms controls:
%listFile = '/data/shimony/shimony2/wunder/wunder_caf_II/ABCD_DCAN_docker/docker_output_BP_filter_and_cleaning_2020_05_17/dconns/N52_termcontrols_subjdconns.lst';
% preemies:
listFile = '/data/shimony/shimony2/wunder/wunder_caf_II/ABCD_DCAN_docker/docker_output_BP_filter_and_cleaning_2020_05_17/dconns/N56_preemies_subjdconns.lst';

% Output directory (passed into your function)
outputdir = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/net_variant_analysis/spatial_corr_maps/output_910_data';

% Read the file
fid = fopen(listFile, 'r');
subDconnLocs = textscan(fid, '%s');
fclose(fid);

subDconnLocs = subDconnLocs{1};  % extract cell array

% Loop through subjects
for i = 1:length(subDconnLocs)

    subDconnLoc = subDconnLocs{i};

    fprintf('Processing %d/%d:\n%s\n', ...
        i, length(subDconnLocs), subDconnLoc);

    try
        createSpatialCorrMap_MR(subDconnLoc, outputdir);
    catch ME
        warning('Failed for %s', subDconnLoc);
        disp(ME.message);
    end

end
