function myRadio(RadioH, EventData)
handles = guidata(RadioH);
otherRadio = handles.radio(handles.radio ~= RadioH);
set(otherRadio, 'Value', 0);
end