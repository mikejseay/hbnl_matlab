function adorn_plots(row_labels,column_labels,xlabel,ylabel,overtitle,subplot_dims)

% plot subplot row/column labels and an "over-title"
ax = findobj(gcf,'type','axes');
pos = cell2mat(get(ax,'position'));
if nargin < 6
    subplot_dims(2) = numel(unique(pos(:,1))); % same X positions
    subplot_dims(1) = numel(unique(pos(:,2)));
end

for row=1:subplot_dims(1)
    subplot(subplot_dims(1),subplot_dims(2),subplot_dims(2)*(row-1)+1);
    rowlabel_locs=determine_labellocs('row');
    hold on; text(rowlabel_locs(1),rowlabel_locs(2),row_labels{row},'FontWeight','bold'); hold off;
end

for column=1:subplot_dims(2)
    subplot(subplot_dims(1),subplot_dims(2),column);
    columnlabel_locs=determine_labellocs('column');
    hold on; text(columnlabel_locs(1),columnlabel_locs(2),column_labels{column},'FontWeight','bold'); hold off;
end

plottitle(overtitle);
plotlabel(xlabel,ylabel)

end