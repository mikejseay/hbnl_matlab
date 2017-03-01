function makescale(pos,vert,horz,horz_text)

%horz line & text
line([pos(1) pos(1)+horz],[pos(2) pos(2)],'Color','k');
text(pos(1)+horz/8,pos(2)-vert/2,horz_text);

%vert line & text
line([pos(1) pos(1)],[pos(2) pos(2)+vert],'Color','k');
text(pos(1)-horz*.75,pos(2)*1.25,num2str(vert));

end