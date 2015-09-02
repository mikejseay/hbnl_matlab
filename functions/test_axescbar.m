%understand axes and colorbars

figure;
a=axes('Position', [0.88 0.05 0.1 0.9], 'Visible', 'off');
set(a,'CLim',[-pi pi]);
c=colorbar('YTick',[-pi,0,pi],'YTickLabel',{'0',[char(177),'pi/2'],[char(177),'pi']});
