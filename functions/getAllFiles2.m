% I used the code mentioned in this great answer and expanded it to support 2 additional parameters which I needed in my case. 
% The parameters are file extensions to filter on and a flag indicating whether to concatenate the full path to the name of the file or not.
% I hope it is clear enough and someone will finds it beneficial.

% Example for running the code:

    % fileList = getAllFiles2('dirName', 'ern_*_avg.h1', 1); % [0=false; 1=true]

function fileList = getAllFiles2(dirName, fileExtension, appendFullPath)
  dirInput=[dirName fileExtension];
  dirData = dir(dirInput);      %# Get the data for the current directory
  dirWithSubFolders = dir(dirName);
  dirIndex = [dirWithSubFolders.isdir];  %# Find the index for directories
  fileList = {dirData.name}';  %'# Get a list of the files
  if ~isempty(fileList)
    if appendFullPath
      fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
    end
  end
  subDirs = {dirWithSubFolders(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles2(nextDir, fileExtension, appendFullPath)];  %#ok<*AGROW> %# Recursively call getAllFiles
  end

end