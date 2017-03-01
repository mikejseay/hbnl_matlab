function session = identify_cntsession_folder(folder)

contents = dir(folder);
n_files = length(contents);
file_ind = 3; % start with 3rd file in dir
session = 'x'; % initialize session identified

% search through each file
while ~ismember(session, {'a','b','c','d','e','f','g','h','i','j','k'}) && ...
        file_ind < n_files
        
    afile = contents(file_ind).name;
    if length(afile) > 6
        session = afile(7); % 7th letter in filename
    end
    file_ind = file_ind + 1;
end

if ~ismember(session, {'a','b','c','d','e','f','g','h','i','j','k'})
    fprintf('Could not identify the session from the folder');
end

end