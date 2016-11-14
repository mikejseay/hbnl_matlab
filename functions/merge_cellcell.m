function out = merge_cellcell( in )

col_len = length( in );
row_len = unique( cellfun(@length, in) );

out = cell(row_len, col_len);

for c = 1:col_len
    
    out(:, c) = in{c};
    
end

end