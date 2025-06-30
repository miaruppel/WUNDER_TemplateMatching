function params = make_fs_masks_params_NU

project_dir = '/projects/b1081/';
singularity_cmd_start = 'module load singularity; ';
afni_sif = [project_dir '/singularity_images/afni_latest.sif'];
templateflow_dir = [project_dir '/Atlases/templateflow/'];
fsl_cmd_start = 'module load fsl/5.0.8; ';
ants_sif = [project_dir '/singularity_images/ants_v2.4.0.sif'];

%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end

end
