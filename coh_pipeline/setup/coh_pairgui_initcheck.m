function startover=coh_pairgui_initcheck(pairs)

n_pairs=size(pairs,1);

geometry = { [3 1] ...
            };

uilist = {  { 'Style', 'text', 'string', ['Your choice resulted in ',num2str(n_pairs),' pairs. To proceed to plot/filter pairs, check the box. Leave unchecked to tweak initial parameters.'], 'horizontalalignment', 'right' }, ...
            { 'Style', 'checkbox', 'Value', 0, 'string', '', 'horizontalalignment', 'right'}, ...
            };

results=inputgui(geometry,uilist,'help coh_pairgui_initcheck;','Check pairs');

startover = ~logical(results{1});

end