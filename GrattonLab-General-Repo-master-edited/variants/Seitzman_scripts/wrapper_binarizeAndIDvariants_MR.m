% Location of spatial correlation maps
% term controls:
inputdir = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/net_variant_analysis/spatial_corr_maps/output_910_data/termcontrols/ABCD_template';
% preemies:
%inputdir = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/net_variant_analysis/spatial_corr_maps/output_910_data/preemies/ABCD_template';

% (Optional) SNR mask 
snrMask = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/net_variant_analysis/bottomBrainMaskFixed.dtseries.nii';

% Output directory
% term controls:
outputdir = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/net_variant_analysis/binarizedIDs/output_910_data/termcontrols/ABCD_template';
% preemies: 
%outputdir = '/data/smyser/smyser4/wunder/wunder_caf_III/TemplateMatching/net_variant_analysis/binarizedIDs/output_910_data/preemies/ABCD_template';

% Get list of spatial corr maps
maps = dir(fullfile(inputdir, '*.dtseries.nii'));
fprintf('Found %d spatial correlation maps\n', length(maps));

% Loop through subjects
for i = 1:length(maps)

    spatialCorrMap = fullfile(maps(i).folder, maps(i).name);

    fprintf('Processing %d/%d:\n%s\n', ...
        i, length(maps), spatialCorrMap);

    try
        if isempty(snrMask)
            binarizeAndIDvariants_MR(spatialCorrMap, [], outputdir);
        else
            binarizeAndIDvariants_MR(spatialCorrMap, snrMask, outputdir);
        end
    catch ME
        warning('FAILED for %s', spatialCorrMap);
        disp(ME.message);
    end
end
