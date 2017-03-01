function pairs=bipartition(fromset,toset)

n_from=length(fromset);
n_to=length(toset);

pairs=zeros(n_from*n_to,2);

dum=0;
for f=1:n_from
    for t=1:n_to

        dum=dum+1;
        
        pairs(dum,:)=[fromset(f),toset(t)];
    end
end

end