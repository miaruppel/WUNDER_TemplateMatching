% Code to cut up tmask to represent various amounts of data in the same subject 
% i.e., cutting a 20 minute run down to 10 minutes and examining how data quality is affected 
% To find minutes: # of frames x 0.8 / 60

inputFile = '/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/recon_docker/docker_output/sub-neo097f/ses-None/files/DCANBOLDProc_v4.0.0/analyses_v2/motion/sub-neo097f_tmask.txt';
outputFile = '/net/10.20.145.47/SMYSER04/smyser4/wunder/wunder_caf_III/recon_docker/docker_output/sub-neo097f/ses-None/files/DCANBOLDProc_v4.0.0/analyses_v2/motion/sub-neo097f_tmask_20mins.txt';
targetCount = 1500; % 20 minutes 

count = 0;
data = {};

file = fopen(inputFile, 'r');
if file == -1
    error('Cannot open file: %s', inputFile);
end 

lineIndex = 1;
while ~feof(file)
    line = fgetl(file);
    if strcmp(line, '1')
        count = count + 1;
    end 
    
     if strcmp(line, '1') && count == targetCount
        fprintf('Count of %d has been reached at line index %d \n', count, lineIndex);
    end 
    
    if count > targetCount
        line = '0';
    end 
    
    data{lineIndex} = line;
    lineIndex = lineIndex + 1;
end 

dlmwrite(outputFile, data, 'delimiter', '\n');
fclose(file);

