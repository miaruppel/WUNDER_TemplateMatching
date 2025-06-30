function Run_Infomap_GrattonLab(Infomap_params_file)
%Run_Infomap(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, [numpools],[structure_indices])
%
% Run infomap on a matrix with a given distance exclusion, at various
% density thresholds, and write the results from all thresholds into a
% single text file named "rawassn.txt". This can take a long time for large
% matrices. It will run up to eight infomaps simultaneously if the Parallel
% Computing Toolbox is installed.
%
% Inputs:
%
% rmat - a correlation matrix to be infomapped. Can be a numeric matrix or
%  a cifti file that will be loaded
% dmatname - a .mat file containing a node-to-node distance matrix
% xdistance - the distance exclusion to apply, in mm (i.e., nodes closer
%  than xdistance are not allowed to be connected)
% thresholdarray - a vector of thresholds to apply to the matrix. Infomap
%  will be run separately for each threshold.
% makebinary - whether the matrix is binarized after thresholding. 1 = make
%  it binary; 0 = leave it weighted.
% outdir - the folder results will be written to. Will be created if it
%  doesn't exist.
% numpools - an OPTIONAL scalar input specifying the number of parallel
%  pools to use. Omit or leave empty ([]) to use the default of 6.
%
% structure_indices - an OPTIONAL vector input of the same length as one
%  dimension of rmat that defines known divisions of interest in the matrix.
%  For example, one might provide a vector defining cortical nodes as "1"
%  and subcortical nodes as "2". If this input is provided, connection
%  values in the matrix will be thresholded separately for each combination
%  of defined structures. This is useful if one wants to e.g. allow
%  subcortical nodes to participate in communities even though they 
%  have very weak connection strengths.
%
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
%EMG 06/25/15
% CG - 07/2020 - edits to work at NU



%LOAD INPUT PARAMETERS
[paramspath,paramsname,paramsextension] = fileparts(Infomap_params_file);
origpath = pwd;
if ~isempty(paramspath)
    cd(paramspath)
end
%Load parameters
params = feval(paramsname);
varnames = fieldnames(params);
for i = 1:length(varnames)
    evalc([varnames{i} ' = params.' varnames{i}]);
end
clear varnames params
cd(origpath)

if ~exist('numpools') || isempty(numpools)
    numpools = 6;
end

%Compute correlation matrix for input to Infomap 
rmat = make_corrmat(datafile, ciftiDir);

% if ischar(rmat) %&& length(rmat)>10 && strcmp(rmat(end-9:end),'.dconn.nii')
%     rmat = ft_read_cifti_mod(rmat);
%     rmat = rmat.data;
% end %elseif ischar(rmat) && length(rmat)>4 && strcmp(rmat(end-3:end),'.mat')
% rmat = single(rmat);

% warning off
fprintf('Applying distance threshold\n')
tic

if ischar(dmatname)
    dmat = smartload(dmatname);
elseif isnumeric(dmatname)
    dmat = dmatname;
    clear dmatname
end

% apply a distance exclusion
if ~isempty(xdistance)
    if numel(dmat) > numel(rmat)
        warning('Distance matrix size exceeds correlation matrix size!\n\tCorr mat: %d by %d and distance mat: %d by %d\n\tAssuming dmat is whole-brain and rmat is cortex-only...\n\tFitting dmat to corr mat...',...
            size(rmat,1),size(rmat,2),size(dmat,1),size(dmat,2));
        dmat = dmat(1:size(rmat,1),1:size(rmat,1));
    elseif numel(rmat) > numel(dmat)
        error('Correlation matrix is larger than distance matrix - fix!')
    end
    if isnumeric(xdistance) && (xdistance>=0)
        rmat(dmat < xdistance) = 0;
    else
        error('xdistance is not >=0 or is not numeric.\n');
    end
end

clear dmat
toc

warning off
numanalyses = length(thresholdarray);
numnodes = size(rmat,1);

outdir = [outdir '/' subject];
if ~exist(outdir)
    mkdir(outdir);
end

dlmwrite([outdir '/thresholds.txt'],thresholdarray,'delimiter',' ')

fprintf('Form matrix and find r thresholds\n');
disp(' ')
%Remove diagonals
for i = 1:numnodes
    rmat(i,i) = 0;
end

tic
if exist('structure_indices','var') && (numel(unique(structure_indices)) > 1)
    ind = matrix_thresholder_structurespecific(rmat,thresholdarray(end),structure_indices);
else
    ind = matrix_thresholder_simple(rmat,thresholdarray(end));
end
toc

if makebinary
    rmat=ceil(rmat);
end

% Save out largest pajek file, i.e. top threshold
pajekfileorig = [outdir '/pajek_col' num2str(length(thresholdarray)) '.net'];
mat2pajek_byindex(rmat,ind,pajekfileorig)
numpossibleedges=(numnodes*(numnodes-1))/2;

clear rmat

% Write out other thresholds
for i = 1:numanalyses
    if i<numanalyses
        pajekfile = [outdir '/pajek_col' num2str(i) '.net'];
        edgesleft=round(thresholdarray(i)*numpossibleedges);
        numuse = edgesleft + numnodes + 2; % Number of edges plus number of nodes plus 2 lines for the headers in the pajek file
        %evalc(['!head -n ' num2str(numuse) ' ' pajekfileorig ' >! ' pajekfile]); %CG edits to make it work
        system(['head -n ' num2str(numuse) ' ' pajekfileorig ' > ' pajekfile]);
    end
end
    


% do analyses at each threshold
%eval(['matlabpool open ' num2str(numpools)]) % CG - turned this off for now
pool=parpool(numpools);  %2/8/22 added this back in
parfor i=1:numanalyses
    fprintf('Thr/box %d, pass %d\n',i);
    
    tic
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    rawclrs = run_infomap_on_pajekfile(pajekfile,100);
    dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
    toc
end

for i = 1:numanalyses
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end



for i = 1:numanalyses
    rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
end
% write the raw assignments as .txt
dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
delete([outdir '/rawassn_col*.txt'])

delete(pool) %close parallel workers
    
end

function rmat = make_corrmat(datafile, ciftiDir)

% LOAD SUBJECT DATALIST
df = readtable(datafile); %reads into a table structure, with datafile top row as labels
numdatas=size(df.sub,1); %number of datasets to analyses (subs X sessions)


fprintf('Loading data and masking with tmask...\n')
catData = []; %stores data for each run with tmask applied

for i = 1:numdatas %For each task in each session
    
    if isa(df.runs(i),'double') && length(df.runs(i))==1
        runnums = table2array(df{i,'runs'});
    else
        runnums = str2num(df.runs{i});
    end
    
    if isa(df.sub(i), 'double')
        sub = num2str(df.sub(i));
    elseif isa(df.sub(i), 'cell')
        sub = df.sub{i};
    end
    
    for j = 1:length(runnums) %For each run in each task

        conf_fstring = sprintf('%s/%s/fmriprep/sub-%s/ses-%d/func/',df.topDir{i},df.confoundsFolder{i},sub,df.sess(i));
        all_fstring = sprintf('sub-%s_ses-%d_task-%s_run-%d',sub,df.sess(i),df.task{i},runnums(j));

        % LOAD DATA        
        ciftifile = sprintf('%s/sub-%s/ses-%d/cifti_timeseries_normalwall/sub-%s_ses-%d_task-%s_run-%d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',ciftiDir,sub,df.sess(i),sub,df.sess(i),df.task{i},runnums(j));
        tmp = ft_read_cifti_mod(ciftifile);
        input_data = tmp.data;
        clear tmp;
        
        % LOAD TMASK
        tmask_fname = [conf_fstring 'FD_outputs/' all_fstring '_desc-tmask_' df.FDtype{i} '.txt'];
        tmask_data = table2array(readtable(tmask_fname));
        tmask_data = logical(tmask_data);
        
        % APPLY tmask to data
        data_masked = input_data(1:59412,tmask_data); %assume input data is grayoordinate X time
        catData = [catData data_masked];
        clear input_data;
    
    end
end
    
% Take correlations among vertices
% a single run takes about 1 min
%disp('Starting dconn correlations...')
fprintf('Running correlations on data size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
tic;
rmat = paircorr_mod(catData');
rmat = single(rmat);
toc;

disp('Finished correlations')

end
