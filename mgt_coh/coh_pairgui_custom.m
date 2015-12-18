function [pairs,hyp_inds,hyp_indlbls]=coh_pairgui_custom

% custom .mat should contain pairs, hyp_inds, and hyp_indlbls

commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a file'');' ...
                    'if filename(1) ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];

geometry = { [3 4 1] 
            };

uilist = {  { 'Style', 'text', 'string', '.mat file containings ''pairs'', ''hyp_inds'', and ''hyp_indlbls'' variable', 'horizontalalignment', 'right' }, ...
              { 'Style', 'edit', 'string', '/export/home/mike/matlab/origin/coords/pairs_FPO92.mat', 'horizontalalignment', 'left', 'tag',  'Pairpath' }, ...
              { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''Pairpath'';' commandload ] }, ...
            };

results=inputgui(geometry,uilist,'help coh_pairgui_custom;','Use custom pairs from a .mat file');

pairfile=results{1};
load(pairfile)

end