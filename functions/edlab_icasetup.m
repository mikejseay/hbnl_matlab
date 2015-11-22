function edlab_icasetup(filename,timeepochvector, varargin)
%% EDLAB_ICASETUP - Imports CNT/EEG data into EEGLAB and sets it up for ICA
%  
%GUI Usage;
%  >> edlab_icasetup;
%
%Non-GUI Useage (useful for script writing): 
%  >> edlab_icasetup(filename,[EpochStart EpochEnd], 'optionalinput1', 
%     'valueofinput1', 'optionalinput2', 'valueofinput2' , ...);
%
%
%Non-GUI Required Inputs: 
%  filename - Neuroscan .cnt or .eeg file
%
%  [EpochStart EpochEnd] - 2 element vector of epoch limits 
%  in seconds.  
%
%
%Non-GUI Optional Inputs: 
%   'Neuropath' (string), path (folder) to locate Neuroscan data. Default is current directory (cd).
%
%   'Newpath' (string), path (folder) to save new files. Default is current directory (cd).
%
%   'Newname' (string), name that you wish to give to the new files, no file extension needed.
%   Default is an acronym from the original Neuroscan file name that
%   recognizes capital A-Z, 1-9 and the underscore character.
%
%   'Map' (string), electrode map file name. Defaults are '64 ch_EEGLab.ced' for .cnt 
%   files and '68 ch_EEG.ced' for .eeg files
%
%   'Mappath' (string), electrode map path (folder). Default is current directory (cd).
%
%   'InclExcl' (integer, 1 or 0). 1 indicates the following electrodes will
%   be included (see directly below); 0 indicates the following electrodes will be
%   excluded. Default is to exclude the listed (or unlisted) electrodes
%   (0). 1 should not be entered unless a list for 'Electrodes' is given,
%   as defaults for 'Electrodes' pertain to excluded electrodes (see
%   directly below).
%
%   'Electrodes' (set of strings), excluded/included electrodes from electrode map. Default for '64
%   ch_EEGLab.ced' is {'HEOG' 'VEOG'} and '68 ch_EEGLab.ced' is {'HEOG' 'VEOG' 'GFP' 'REF'}
%
%   'Resample' (integer), rate at which to resample the data. Default is set at 250 Hz.
%
%   'Behavioralname' (string), file name of behavioral data for behavioral merge.
%   If no name is given, behavioral merge will not be performed
%
%   'Behavioralpath' (string), path (folder) to locate behavioral data for behavioral merge. 
%   Default is current directory (cd)
%
%   'NoRTcode' (real number), when a behavioral merge is formed, code for no reaction time.
%   Default is to leave the code blank ( [] ).
%
%
%Usage Examples:
%  >> edlab_icasetup
%       Opens the GUI window.
%
%  >> edlab_icasetup('Michael Nunez.cnt',[-.5 1])
%       Sets up EEGLAB data from the Neuroscan file 'Michael Nunez.cnt' with 
%       epochs starting at -.5 seconds before stimulus and ending 1.5 seconds 
%       after stimulus
%
%  >> edlab_icasetup( 'Michael Nunez.cnt' , [-.5 1] , 
%  'Neuropath' , 'C:\Users\Dr. Jekyll\Study 1\Neuroscan Data\' , 
%  'Newpath' , 'C:\Users\Dr. Jekyll\Study 1\EEGLAB Data\',
%  'Newname' , 'Subject 1'
%  'Map' , '64 ch_EEGLab.ced',
%  'Mappath' , 'C:\EEGLAB\Electrode locations\' ,
%  'InclExcl' , 0 ,
%  'Electrodes' , {'HEOG' 'VEOG' 'PO3'} , 
%  'Resample', 500 , 
%  'Behavioralname' , 'MNBehavior.dat' , 
%  'Behavioralpath' , 'C:\Users\Mr. Hyde\Study 1\Neuroscan Data\' , 
%  'NoRTcode' , 2000 )
%       This is an example of how to use the non-gui version of the function in
%       its extended form.
%
%Notes:
%  - The defaults follow standard procedure of Dr. Golob's Lab as of 9/14/11
%  - For non-GUI, it is easiest to put all of the files that this function uses in a
%  folder and change MATLAB's current directory to that folder; the default
%  paths are all set to 'cd', MATLAB's current directory.
%  - For non-GUI, a '.ced' file is always needed in your current directory, or one must
%  specifiy a 'Map' and/or 'Mappath'
%  - For non-GUI, there does not have to be any particular order to optional inputs

%% Record of revisions:
%   Date           Programmers               Discription of change
%   ====        =================            =====================
%  9/6/11         Michael Nunez       Original code using EEGLAB functions
%                                           with use of John DiLeo 
%                                         and Steuart Turner scripts.
%
%  9/14/11        Michael Nunez     Added GUI option using an EEGLAB function.  
%                                   Added the option to Include or Exclude
%                                             electrodes.

%% GUI
if nargin < 2
    
    % popup window parameters
    % -----------------------

commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];

dirload = [ '[dirpath] = uigetdir([],''Location to save new files.'');' ...
                    'if dirpath(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [dirpath]);' ...
                    'end;' ...
                    'clear dirpath tagtest;' ];

geometry    = { [2 1 1.5] [2 1 1.5] [1] [2 1 1.5] [3.33 1.66 1 1.5] [2 1 1.5] [2 1 1.5] [2 1 1.5] [2 1 1.5] ...
                [1] [2 1 1.5] [2 1 1.5]};
                
uilist = {    { 'Style', 'text', 'string', 'Open Neuroscan .cnt or .eeg file.', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'Neuroscanfile' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Neuroscanfile'';' commandload ] }, ...
              { 'Style', 'text', 'string', 'Open electrode map.', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'Mapfile' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Mapfile'';' commandload ] }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', 'Checked -> Include    Unchecked -> Exclude', 'fontweight', 'bold'  }, ...
              { 'Style', 'text', 'string', 'Channels to exclude or include:' }, ...
              { 'Style', 'edit', 'string', '', 'tag', 'chans' }, ...
              {}, { 'Style', 'checkbox', 'string', '' , 'horizontalalignment', 'center'}, ...
              { 'Style', 'text', 'string', 'Sampling rate:', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '250' }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', 'Epoch limits [start end], in seconds:', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '[-.5 1]' }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', 'Open behavioral data (if applicable).', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'Datfile' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Datfile'';' commandload ] }, ...
              { 'Style', 'text', 'string', 'Integer code for no reaction time (applicable only for behavioral merge).', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', ''  }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', ''  }, ...
              { 'Style', 'text', 'string', 'Name of new files (no path or file extension included).', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '' }, ...
              { 'Style', 'text', 'string', '' }, ...
              { 'Style', 'text', 'string', 'Folder to place new files:', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'Newdirpath' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Newdirpath'';' dirload ] } };
          
results = inputgui(geometry, uilist, 'help edlab_icasetup;', 'Set Up Data for ICA Analysis - edlab_icasetup');

%Makes to see all required inputs were given

if isempty(results{1}) == 1
    error('A .cnt or .eeg Neuroscan file must be given.');
end

if isempty(results{2}) == 1
    error('An electrode map must be given.');
end

if isempty(results{5}) == 1
    error('A resampling rate must be given.');
end

if isempty(results{6}) == 1
    error('Epoch limits must be given.');
end

    
 %Define variables from GUI
 
   %Separates Neuroscan file name and Neuroscan file path
filenamepath = regexp(results{1}, '\\', 'split');
filename = filenamepath{length(filenamepath)};
opt.Neuropath = 'C:';
for n = 2:(length(filenamepath) - 1)
    opt.Neuropath = [opt.Neuropath, '\', filenamepath{n}];
end
opt.Neuropath = [opt.Neuropath, '\'];
   
    %Separates Electrode Map file name from path
mapnamepath = regexp(results{2}, '\\', 'split');
opt.Map = mapnamepath{length(mapnamepath)};
opt.Mappath = 'C:';
for n = 2:(length(mapnamepath) - 1)
    opt.Mappath = [opt.Mappath, '\', mapnamepath{n}];
end
opt.Mappath = [opt.Mappath, '\'];

    %Finds accepted or rejected electrodes
opt.Electrodes = regexp(results{3}, '\s', 'split');

    %Accept or reject electrodes?
opt.InclExcl = results{4};

    %Finds resample rate
opt.Resample = str2num(results{5});

    %Finds epoch boundaries
timeepochvector = str2num(results{6});

    %Separates Behavioral file name from path, only if applicable
opt.Behavioralpath = [];
if isempty(results{7}) == 1
    opt.Behavioralname = 'nnvphvgueanrnmh';
else
    behavioralnamepath = regexp(results{7}, '\\', 'split');
    opt.Behavioralname = behavioralnamepath{length(behavioralnamepath)};
    opt.Behavioralpath = 'C:';
    for n = 2:(length(behavioralnamepath) - 1)
       opt.Behavioralpath = [opt.Behavioralpath, '\', behavioralnamepath{n}];
     end
opt.Behavioralpath = [opt.Behavioralpath, '\'];
end

    %Finds code for no RT
opt.NoRTcode = str2num(results{8}); 
    
    %Finds name for new files
opt.Newname = results{9};

    %Defines path for new files
if isempty(results{10}) == 1
    opt.Newpath = [cd, '\'];
else
    opt.Newpath = [results{10}, '\'];
end

    
    
end
%% Preliminary Code

%This will be used to detect whether Neuroscan file is .cnt or .eeg
CNTorEEG = regexp(filename, '\.', 'split');

%Returns an acronym from the CNT or EEG file name which will be used as the
%subject's initials
capitals = regexp(filename, '[A-Z0-9_]', 'match');
acronym = [];
for n = 1:length(capitals)
    acronym = [acronym , capitals{n}];
end

%Only for GUI, set the new file name to the acronym if no name was given
if nargin < 2
    if isempty(results{9}) == 1
    opt.Newname = acronym;
    end
end
%% Defining Optional Inputs For Non-GUI

if nargin >= 2

switch CNTorEEG{2}
    case 'cnt' %Do the following steps for .cnt Neuroscan files
        
    %Sets electrode map default for .cnt files
    MapDefault = '64 ch_EEGLab.ced';
    %Sets excluded electrodes default for .cnt files
    ElectrodesDefault = {'HEOG' 'VEOG'};

    case 'eeg' %Do the following steps for .eeg Neuroscan files

    %Sets electrode map default for .eeg files
    MapDefault = '68 ch_EEGLab.ced';
    %Sets excluded electrodes default for .eeg files
    ElectrodesDefault = {'HEOG' 'VEOG' 'GFP' 'REF'};
    
    otherwise
    error('Neuroscan file must be either in .cnt or .eeg format. If it is, please include the file extension in the string input.');

end


opt = finputcheck( varargin, { 'Neuropath'       'string'  [] [cd,'\'];
                               'Newpath'         'string'  [] [cd,'\'];
                               'Newname'         'string'  [] acronym;
                               'Map'             'string'  [] MapDefault;
                               'Mappath'         'string'  [] [cd,'\'];
                               'InclExcl'        'integer' [] 0;
                               'Electrodes'      'cell'    [] ElectrodesDefault;
                               'Resample'        'integer' [] 250;
                               'Behavioralname'  'string'  [] 'nnvphvgueanrnmh';
                               'Behavioralpath'  'string'  [] [cd,'\'];
                               'NoRTcode'        'real'    [] NaN } );
                               
if isstr(opt), error(opt); end;

if isequal(opt.Electrodes , ElectrodesDefault)
switch opt.Map
    case '64 ch_EEGLab.ced'
    
    %Sets excluded electrodes default for '64 ch_EEGLab.ced' channel location file
    opt.Electrodes = {'HEOG' 'VEOG'};
    
    case '68 ch_EEGLab.ced'
    
     %Sets excluded electrodes default for '68 ch_EEGLab.ced' channel location file
    opt.Electrodes = {'HEOG' 'VEOG' 'GFP' 'REF'};
    
end
end

if opt.InclExcl ~= 0 & opt.InclExcl ~= 1
    error('Only an integer (non-string) value of 0 or 1 is appropriate for ''InclExcl');
end

end
%% EEGLAB Code

%(Re)opens EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

if all(CNTorEEG{2} == 'cnt') %Do the following steps for .cnt Neuroscan files
    
    %Imports .cnt dataset, sets data format to 32 bit
    EEG = pop_loadcnt([opt.Neuropath, filename] , 'dataformat', 'int32');

elseif all(CNTorEEG{2} == 'eeg') %Do the following steps for .eeg Neuroscan files

    %Imports .eeg dataset, sets data format to 32 bit
    EEG = pop_loadeeg(filename, opt.Neuropath, [],[],[],[], 'int32');
    
else
    
    error('Neuroscan file must be either in .cnt or .eeg format.');

end

%Save dataset
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew', [opt.Newpath , opt.Newname ,'.set'] ,'gui','off');
EEG = pop_loadset('filename', [opt.Newname,'.set'] ,'filepath', opt.Newpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%Get channel locations from channel location file.
EEG = pop_chanedit(EEG, 'lookup', [opt.Mappath,opt.Map]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

if opt.InclExcl == 0
%Exclude the list of electrodes
EEG = pop_select( EEG,'nochannel', opt.Electrodes);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
EEG = eeg_checkset( EEG );
elseif opt.InclExcl == 1
%Include the list of electrodes
EEG = pop_select( EEG, 'channel', opt.Electrodes);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
EEG = eeg_checkset( EEG );
end

%Bandpass filter between 1 and 50 hz
EEG = pop_iirfilt( EEG, 1, 0, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
EEG = pop_iirfilt( EEG, 0, 50, [], [0]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
EEG = eeg_checkset( EEG );

%Resample the dataset
EEG = pop_resample( EEG, opt.Resample);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );

%Reference the data using the average reference
EEG = pop_reref( EEG, []);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
EEG = eeg_checkset( EEG );

%Remove channel baseline across entire time period
EEG = pop_rmbase( EEG, [EEG.xmin*1000 EEG.xmax*1000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%Save and name unepoched dataset
EEG = pop_saveset( EEG, 'filename',[opt.Newname,'_unepoched.set'],'filepath',opt.Newpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

%Extract data epochs
EEG = pop_loadset('filename', [opt.Newname,'_unepoched.set'] ,'filepath', opt.Newpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
EEG = pop_epoch( EEG, {  }, timeepochvector, 'newname', 'Epoched Data', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
%EEG = pop_rmbase( EEG, [], [] );
EEG = pop_rmbase( EEG, [EEG.xmin*1000 EEG.xmax*1000]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%Save and name epoched dataset
EEG = pop_saveset( EEG, 'filename',[opt.Newname,'_epoched.set'],'filepath',opt.Newpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

%Merge behavioral data if applicable

switch opt.Behavioralname
    case 'nnvphvgueanrnmh' %Case when no Behavioral data file name is given, behavioral merge is not performed. 'nnvphvgueanrnmh' is a meaningless string of characters
    otherwise
EEG = pop_loaddat( EEG , [opt.Behavioralpath,opt.Behavioralname], opt.NoRTcode);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%Save and name epoched dataset with behavioral data if applicable
EEG = pop_saveset( EEG, 'filename',[opt.Newname,'_epoched_behavioral.set'],'filepath',opt.Newpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

end

disp('Done.');
