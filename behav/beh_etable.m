function behdata = beh_etable(etable, expstruct)
% interpret etable

stimtypes   = expstruct.ttl_codes;
n_types     = length(stimtypes);

%initialize vars
[ acc, accwithlate, medianrt, medianrtwithlate] = deal(nan(n_types, 1));

% for each stimulus type (ttl code)

for stimtype = 1:n_types
    
    % get stims
    tmp_stimevs = etable.type == stimtypes( stimtype );
    
    % if any response was required
    if ~ismember( expstruct.correct_resps(stimtype), [0 -1] )
        % resps should be immediately following
        % correct_resp = expstruct.correct_resps(stimtype);
        tmp_respevs = circshift(tmp_stimevs, 1) & ...
            etable.resp_code ~= 0;

        % correct (excludes late responses)
        tmp_correct = tmp_respevs & etable.correct & ~etable.late;
        % correct (includes late responses)
        tmp_correct_withlate = tmp_respevs & etable.correct;
        
        acc(stimtype)               = sum(tmp_correct) ./ sum(tmp_stimevs);
        accwithlate(stimtype)      = sum(tmp_correct_withlate) ./ sum(tmp_stimevs);
        medianrt(stimtype)          = nanmedian(etable.rt(tmp_correct));
        medianrtwithlate(stimtype) = nanmedian(etable.rt(tmp_correct_withlate));
    
    elseif expstruct.correct_resps(stimtype) == 0 % if absence of response was required
        tmp_correct = etable.correct(tmp_stimevs);
        
        acc(stimtype)               = sum(tmp_correct) ./ sum(tmp_stimevs);
        
    end % if there is no accuracy/RT information that is meaningful
        % leave everything as NaNs
end

name = expstruct.name;
ttl_descriptors = expstruct.ttl_descriptors;

behdata = v2struct(name, ttl_descriptors, acc, accwithlate, medianrt, ...
    medianrtwithlate, etable);

end

%{
meanrt(stimtype)            = nanmean(etable.rt(tmp_correct));
meanrt_withlate(stimtype)   = nanmean(etable.rt(tmp_correct_withlate));
stdrt(stimtype)             = nanstd(etable.rt(tmp_correct));
stdrt_withlate(stimtype)    = nanstd(etable.rt(tmp_correct_withlate));
%}