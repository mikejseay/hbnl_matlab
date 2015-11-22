function [pairs,hyp_inds,hyp_indlbls]=coh_choosepairs(type,coords_file)

init_startover=true;
while init_startover %initial pair set

switch type
    case 'determine'
        
        pairs=coh_pairgui_determine(coords_file);
        hyp_inds=[];
        hyp_indlbls={};        
        
    case 'seed'
        
        [pairs, hyp_inds, hyp_indlbls]=coh_pairgui_seed;
        
    case 'region'
        
        [pairs, hyp_inds, hyp_indlbls]=coh_pairgui_region(coords_file);
        
    case 'custom'
        
        [pairs, hyp_inds, hyp_indlbls]=coh_pairgui_custom;
        
end
init_startover=coh_pairgui_initcheck(pairs);
end
%initial pair set chosen


second_startover=true;
while second_startover %filter and plot

    [second_startover,plot_bool,filter_bool,restore_bool]=coh_pairgui_check(pairs);
    chan_locs=pass_chan_locs(coords_file);
    
    if plot_bool
        plot_pairs(pairs, chan_locs, hyp_inds);
    end
    
    if filter_bool && ~restore_bool && second_startover
        savepairs=pairs;
        saveinds=hyp_inds;
        [pairs,hyp_inds]=coh_pairgui_filter(pairs,hyp_inds,chan_locs);
    end
    
    if restore_bool && second_startover
        pairs=savepairs;
        hyp_inds=saveinds;
    end

end

end