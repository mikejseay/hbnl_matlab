mgt_coh_list_small=cell(s,1);
for chosen_s=1:32
if s_inds(chosen_s)
mgt_coh_list_small{chosen_s}={};
else
mgt_coh_list_small{chosen_s}=mgt_coh_list{chosen_s};
end
end