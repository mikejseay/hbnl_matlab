function behdata = beh_calc_dat( datpath )
% given a dat file, extract behavioral info using expstruct

% initialize behdata struct
behdata = [];

% build expstruct
expstruct = build_expstruct;

% get index of experiment to target
[~, filename] = fileparts(datpath);
exp_name = filename(1:3);
exp_ind = strcmpi(exp_name, {expstruct.name});
   
try
    % extract event table from CNT
    etable = dat_getbehav( datpath, expstruct(exp_ind) );

    % interpret behavioral data
    behdata = beh_etable( etable, expstruct(exp_ind) );

catch err
    % some reading error most likely (corrupt .mat)
    fprintf('err: %s', err.message);
end

end