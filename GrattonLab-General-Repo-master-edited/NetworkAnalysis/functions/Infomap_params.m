function params = Infomap_params
clear all;


%%% set up data to make correlation matrix
subject = 'INET026';  %will create subject folder within outdir
datafile = '/home/amd3565/INET_Datalists/INET026_Datalist_fmriprep-20.2.0_partial.xlsx';
ciftiDir = '/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI_20.2.0';
outdir = '/projects/b1081/member_directories/adworetsky/dconns2';
dmatname = '/projects/b1081/Scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_large.mat';

%%% define other relevant parameters
xdistance = 30;   %typical is 30 for surface data, 20 for volume data
thresholdarray = [0.003 0.004 0.005:0.005:0.05];
makebinary = 0;
numpools = 8;
structure_indices = [];
cortexOnly = 1; %assumes 59412-vertex (fsLR_32k) surface

if exist(outdir)
    mkdir([outdir '/' subject])
elseif ~exist(outdir,'dir')
    mkdir(outdir)
    mkdir([outdir '/' subject])
end

save([outdir '/' subject '/params_' datestr(now,'mmddyy-HHMMSS') '.mat'], '-v7.3');


%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end
