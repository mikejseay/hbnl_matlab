function linkcaxes(h,sp_inds)
%link c axes to min/max of data given a fig handle and subplot indices

axesObjs = get(h, 'Children');

dataObjs = get(axesObjs, 'Children');

dataObjs = flipud(dataObjs);

dataObjs4 = dataObjs{4};

isprop(dataObjs4, 'ZData')

dataObjs41 = dataObjs4(1);
dataObjs42 = dataObjs4(2);
dataObjs42zdata = get( dataObjs42, 'ZData' );

end