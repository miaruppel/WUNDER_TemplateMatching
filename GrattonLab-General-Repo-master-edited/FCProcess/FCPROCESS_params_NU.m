function params = FCPROCESS_params_NU

project_dir = '/projects/b1081/';
singularity_cmd_start = 'module load singularity; ';
afni_sif = [project_dir '/singularity_images/afni_latest.sif'];
templateflow_dir = [project_dir '/Atlases/templateflow/'];

%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end

end