function out_matlist=match_uniqids(in_matlist,match2_matlist)
%outputs a version of in_matlist that is index-matched to match2_matlist

n_in=length(in_matlist);
n_out=length(match2_matlist);

in_unids=cell(n_in,1);
out_unids=cell(n_out,1);

%extract unique session IDs and make tables for join operation
for file=1:n_in
    [~,fullname,~]=fileparts(in_matlist{file});
    in_unids{file}=fullname(7:17);
end
for file=1:n_out
    [~,fullname,~]=fileparts(match2_matlist{file});
    out_unids{file}=fullname(7:17);
end

[~,indmat,~]=intersect(in_unids,out_unids);

out_matlist={in_matlist{indmat}}';

end