% very easy way to do an x or y crossing

plot(get(gca,'xlim'),[0 0],'k')
plot([6 6],get(gca,'ylim'),'k:')