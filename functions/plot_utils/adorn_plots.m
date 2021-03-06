function adorn_plots(row_labels,column_labels,xlabel,ylabel,overtitle,subplot_dims,exclude_lastcol)

if nargin < 7
    exclude_lastcol=false;
end

% plot subplot row/column labels and an "over-title"
ax = findobj(gcf,'type','axes');
pos = cell2mat(get(ax,'position'));
if nargin < 6
    subplot_dims(2) = numel(unique(pos(:,1))); % same X positions
    subplot_dims(1) = numel(unique(pos(:,2)));
end

if ~isempty(row_labels)
for row=1:subplot_dims(1)
    subplot(subplot_dims(1),subplot_dims(2),subplot_dims(2)*(row-1)+1);
    rowlabel_locs=determine_labellocs('row');
    hold on; text(rowlabel_locs(1),rowlabel_locs(2),row_labels{row},'FontWeight','bold'); hold off;
end
end

if true
    if false %if we're probably plotting time windows
        n_cols = subplot_dims(2) - 1; %skip the last (we get it)
    else
        n_cols = subplot_dims(2);
    end
else
    n_cols = subplot_dims(2);
end

if ~isempty(column_labels)
    if exclude_lastcol
        subtractor=1;
    else
        subtractor=0;
    end
for column=1:n_cols-subtractor
    %subplot(subplot_dims(1),subplot_dims(2),(subplot_dims(1)-1)*subplot_dims(2)+column);
    subplot(subplot_dims(1),subplot_dims(2),column);
    %if column==n_cols
    %    colorbar;
    %end
    columnlabel_locs=determine_labellocs('column');
    hold on; text(columnlabel_locs(1),columnlabel_locs(2),column_labels{column},'FontWeight','bold'); hold off;
end
end

if iscell(overtitle) %quick fix
else
plottitle(overtitle,2);
end
plotlabel(xlabel,ylabel)

end