function opt = beh_createopt
% beh_createopt - GUI for creation of options structure for batch analysis
% of EEG data files in the HBNL lab at SUNY Downstate Medical Center.

%  GUI
if nargin < 2
    
    % popup window parameters
    % -----------------------

commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a file'');' ...
                    'if filename(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];

dirload = [ '[dirpath] = uigetdir([],''Select a folder.'');' ...
                    'if dirpath(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [dirpath]);' ...
                    'end;' ...
                    'clear dirpath tagtest;' ];

geometry    = { [2 1 1.5] ... %[2 1 1 1 1]... %tf specs
                [1] [2 1 1.5] [2 1 1.5] [2 1 1.5]}; %outputs
                
uilist = {    { 'Style', 'text', 'string', 'Text file containing list of full data file path(s)', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '/', 'horizontalalignment', 'left', 'tag',  'Pathfile' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Pathfile'';' commandload ] }, ...
              ...
              %{
              { 'Style', 'text', 'string', 'Measures to include in output .mats' }, ...
              { 'Style', 'checkbox', 'Value', 1 'string', 'Good Trial Series', 'horizontalalignment', 'right'}, ...
              { 'Style', 'checkbox', 'Value', 1 'string', 'ERO (total power)', 'horizontalalignment', 'right'}, ...
              { 'Style', 'checkbox', 'Value', 1 'string', 'ITC', 'horizontalalignment', 'right'}, ...
              { 'Style', 'checkbox', 'Value', 0 'string', 'Coherence', 'horizontalalignment', 'right'}, ...
              %}
              ...
              ... %
              ...
              { 'Style', 'text', 'string', ''  }, ...
              ...
              { 'Style', 'text', 'string', 'Suffix for new .mat files', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '_' }, ...
              { 'Style', 'text', 'string', '' }, ...
              ...
              { 'Style', 'text', 'string', 'Folder to place output .mats:', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '/active_projects/mike/', 'horizontalalignment', 'left', 'tag',  'Outpath' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Outpath'';' dirload ] } ...
              ...
              { 'Style', 'text', 'string', 'Folder to place options file, batch scripts, and logs:', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '/export/home/mike/matlab/batch/', 'horizontalalignment', 'left', 'tag',  'Optpath' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Optpath'';' dirload ] }};
          
results = inputgui(geometry, uilist, 'help beh_createopt;', 'Set Up Option for Data Processing - beh_createopt')';

% Verify inputs
if any(cellfun(@isempty, results([1 end-2:end]) ))
    error('A key input was omitted.');
end

% Define variables from GUI
infile_list         = results{1};
[~,list_filename,~] = fileparts(results{1});

% define expstruct
expstruct           = build_expstruct;

% define exps to examine
exps                = {'vp3','aod','gng','cpt','ern','ant','stp'};

%measure_names   = {'erptrial', 'wave_totpow', 'wave_evknorm', 'coh'};
%measures        = measure_names( find( [results{20:23}] ) );

outsuffix = results{end-2};

outpath = results{end-1};
if ~exist(outpath, 'dir') %make if DNE yet
    mkdir(outpath);
    fprintf('Creating output directory at %s\n', outpath);
end
batchpath = results{end};
if ~exist(batchpath, 'dir')
    mkdir(batchpath);
    fprintf('Creating directory for options, batch, and log files at %s\n', batchpath);
end

%determine the full path where the opt will be saved
optfilename = [list_filename, outsuffix, '_opt.mat'];
optpath     = fullfile(batchpath, optfilename);

%create batch scripts with 100 files per process, or 4 divided equally

%fix the name of the batch_id so that it does not contain illegal chars
batch_id = [list_filename, outsuffix];
if ~isvarname(batch_id)
    batch_id = matlab.lang.makeValidName(batch_id);
end

batchscripts = beh_makebatches(batch_id, infile_list, batchpath, optpath);
n_batches = length(batchscripts);

logpath = cell(n_batches,1);
for batch = 1:n_batches
    logfilename = [list_filename, outsuffix, '_b', num2str(batch), '.log'];
    logpath{batch} = fullfile(batchpath, logfilename);
end

  
opt=v2struct(infile_list, expstruct, exps, ...
    outsuffix, outpath, batchpath, optpath, ...
    batch_id, batchscripts, logpath);

if exist(optpath, 'file')
    overwrite=input('Options file being overwritten. Continue? (y/n) ', 's');
    if strcmpi(overwrite, 'y')
        save(optpath, 'opt');
    else
        fprintf('Exiting...\n');
        return 
    end
else
    save(optpath, 'opt');
end

% make super batches
coh_makesuperbatches(opt,4,2014,'matlab');

disp('Done.');

end