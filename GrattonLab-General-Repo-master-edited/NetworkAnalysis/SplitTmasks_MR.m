% Quickly splitting tmask into individual runs 
% Uses FD trace files to find out accurate run length numbers
% (DICOM.studies file is not always correct!) 

% hard coded subject list 
subjectIDs = {'sub-neo504f','sub-neo539f','sub-caf044f','sub-caf049f','sub-caf059f', ...
    'sub-neo305f','sub-neo510f', 'sub-caf067f','sub-neo299f','sub-neo320f', ...
    'sub-neo314f','sub-neo519f','sub-neo535f','sub-neo508f'};

% there are multiple base directories for these subject IDs (smyser3 and smyser4)
%baseDir = '/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/recon_docker/docker_output';
dirOptions = {'/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/recon_docker/docker_output', ...
    '/net/10.20.145.47/SMYSER03/smyser3/wunder/wunder_caf_III/recon_docker/docker_output'};

% loop through subject IDs
for i = 1:length(subjectIDs)

    subjectID = subjectIDs{i};
    
    % flag to track if files were found
    filesFound = false;
    
    % loop through the directory options
    for k = 1:length(dirOptions)
        
        % relevant directory based on the subject ID and current directory option
        relDir = fullfile(dirOptions{k}, subjectID, 'ses-None', 'files', 'summary_DCANBOLDProc_v4.0.0');
        
        % find files matching the pattern 'FD_task-rest0?.txt'
        filePattern = fullfile(relDir, 'FD_task-rest0*.txt');
        matchingFiles = dir(filePattern);
        
        % check if there are any matching files in this directory
        if ~isempty(matchingFiles)
            filesFound = true;
            fprintf('Found %d FD files for subject %s in directory %s:\n', length(matchingFiles), subjectID, dirOptions{k});
            
            % want to split tmask by run 
            tmaskDir = fullfile(dirOptions{k}, subjectID, 'ses-None', 'files', 'DCANBOLDProc_v4.0.0', 'analyses_v2', 'motion');
            fullTmask = fullfile(tmaskDir, [subjectID, '_tmask_0.1FDThreshold.txt']);

            % load in tmask txt file
            data_tmask = load(fullTmask);
            
            % initialize index for tracking the start position in data_tmask
            startIdx = 1;
            
            % loop through the FD trace files
            for j = 1:length(matchingFiles)
                % full file name
                fullFileName = fullfile(relDir, matchingFiles(j).name);
                
                % open and read file for word count
                FDfile = fopen(fullFileName, 'r');

                % initialize line counter
                count = 0;

                % read file line by line
                while ~feof(FDfile)
                    fgets(FDfile);  % each line
                    count = count + 1;  
                end

                % close file
                fclose(FDfile);

                % each FD trace seemingly substracts -1 to the full run
                % add back this 1
                fullRun = count + 1;
                
                % display the line count
                fprintf('  %s has %d lines.\n', matchingFiles(j).name, fullRun);
                
                % compute the end index for this run 
                endIdx = startIdx + fullRun - 1;
                
                % extract the relevant portion from data_tmask
                runTmask = data_tmask(startIdx:endIdx);
                
                % generate the output file name for the split tmask
                outputTmaskFile = fullfile(tmaskDir, [subjectID, '_task-rest_run-', num2str(j), '-tmask_0.1fFD.txt']);
                
                % save the split tmask data to a new file
                dlmwrite(outputTmaskFile, runTmask, 'delimiter', '\n');
                
                % display info
                fprintf('Created file: %s with %d lines (indices %d to %d).\n', ...
                    outputTmaskFile, length(runTmask), startIdx, endIdx);
                
                % update the start index for the next run
                startIdx = endIdx + 1;
            end
            
            % break out of the directory loop if files were found
            break;
        end
    end
    
    % if no files found in any directory, print a message
    if ~filesFound
        fprintf('No matching FD files found for subject %s in any directory.\n', subjectID);
    end
    
end
    
  