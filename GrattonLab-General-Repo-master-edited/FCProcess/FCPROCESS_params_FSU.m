function params = FCPROCESS_params_FSU

project_dir = '/gpfs/research/grattonlab/';
singularity_cmd_start = '';
afni_sif = [project_dir '/singularity_images/afni_make_build_latest.sif'];
templateflow_dir = [project_dir '/Atlases/templateflow/'];;

%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end

end