function mat = cat_cellcontents(c)

%remove empty cells
%check to make sure mats are same size

[cr, cc] = size(c);

% cat along first dimension (rows)
for col = 1:cc
    c2{col} = cat(1,c{:,col});
end

for row = 1:cr
    mat = cat(2,c2{:});
end

end