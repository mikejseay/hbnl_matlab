%clf
figure(1);
clf

%create panel
p = panel();

%pack config
p.pack({4/15 6/15 5/15});
p(1).pack('h',4);
p(2).pack('h',2);
p(3).pack('h',2);

%check layout
%p.select('all');

%margins
p.de.margin = 3;
p.margintop=15;
p(2).margintop=15;
p(3).margin=[20 15 20 15];
p(3).de.marginright=20;
%p.identify()

%%

%create references
p_scat=p(3,1);
p_hist=p(3,2);

p_scat.select();
%title('Scatter Plot')
xlabel('Phase Diff. (rads)')
ylabel('ITC')
p_hist.select();
%title('Histogram')
xlabel('ITC')


p(2, 1).select();
title('ERP/TF Plot 1')
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
p(2, 2).select();
title('ERP/TF Plot 2')
xlabel('Time (ms)')
%ylabel('Frequency (Hz)')

p(1,1).select();
title('Topoplot 1')
p(1,2).select();
title('Wind Rose 1')

p(1,3).select();
title('Topoplot 2')
p(1,4).select();
title('Wind Rose 2')