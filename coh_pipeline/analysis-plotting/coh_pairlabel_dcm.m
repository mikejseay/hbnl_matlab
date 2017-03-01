function txt = coh_pairlabel_dcm(empt,event_obj)
% Customizes text of data tips

global lookup scl

pos = get(event_obj,'Position');
ind = find(all(bsxfun(@eq,lookup,pos),2));
txt = scl.p_label{ind};