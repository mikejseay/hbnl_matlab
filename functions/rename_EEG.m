function EEG = rename_EEG( EEG, map)
% rename events in an EEG structure

for ev=1:length(EEG.event)
    try
        if isnumeric(EEG.event(ev).type)
            EEG.event(ev).type = map(EEG.event(ev).type);
        elseif ischar(EEG.event(ev).type)
            if strcmpi(EEG.event(ev).type,'rt') || strcmpi(EEG.event(ev).type,'boundary')
                continue
            end
            EEG.event(ev).type = map( str2double(EEG.event(ev).type) );
        end
    catch
        EEG.event(ev).type = 'bad_code';
    end
end

end