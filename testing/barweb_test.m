% understand barweb

a=50+20*randn(100,1);
b=30+10*randn(100,1);
c=60+15*randn(100,1);
d=40+5*randn(100,1);

barvalues=[mean(a),mean(b);mean(c),mean(d)];
errors=[std(a),std(b);std(c),std(d)];

%BASICALLY WHAT THIS MEANS IS THAT BARWEB EXPECTS GROUPS TO BE ROWS
%AND CONDITIONS TO BE COLUMNS

width=[];
groupnames={'CTL';'ALC'};
bw_title='Data';
bw_xlabel=[];
bw_ylabel='Size';
bw_colormap=bone; %pmkmp(length(barvalues),'cubicl');
gridstatus='y';
bw_legend={'Cond1';'Cond2'};
error_sides=1;
legend_type='plot'; %'plot'

figure;
subplot(2,2,1)
bar_h=barweb(barvalues,errors,width,groupnames,bw_title,bw_xlabel,bw_ylabel,...
    bw_colormap,gridstatus,bw_legend,error_sides,legend_type);

clear a b