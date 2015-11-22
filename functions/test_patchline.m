% EXAMPLES:
%% Example 1:
%
figure;
 n = 10;
 xs = rand(n,1);
 ys = rand(n,1);
 zs = rand(n,1)*3;
 plot3(xs,ys,zs,'r.')
 xlabel('x');ylabel('y');zlabel('z');
 p  = patchline(xs,ys,zs,'linestyle','--','edgecolor','g',...
     'linewidth',3,'edgealpha',0.2);

%% Example 2: (Note "hold on" not necessary here!)

 t = 0:pi/64:4*pi;
 p(1) = patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
 p(2) = patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
 l = legend('sine(t)','cosine(t)');
 tmp = sort(findobj(l,'type','patch'));
 for ii = 1:numel(tmp)
     set(tmp(ii),'facecolor',get(p(ii),'edgecolor'),'facealpha',get(p(ii),'edgealpha'),'edgecolor','none')
 end
%
%% Example 3 (requires Image Processing Toolbox):
%%%   (NOTE that this is NOT the same as showing a transparent image on 
%%%         of the existing image. (That functionality is
%%%         available using showMaskAsOverlay or imoverlay).
%%%         Instead, patchline plots transparent lines over
%%%         the image.)
%
 img = imread('rice.png');
 imshow(img)
 img = imtophat(img,strel('disk',15));
 grains = im2bw(img,graythresh(img));
 grains = bwareaopen(grains,10);
 edges = edge(grains,'canny');
 boundaries = bwboundaries(edges,'noholes');
 cmap = jet(numel(boundaries));
 ind = randperm(numel(boundaries));
 for ii = 1:numel(boundaries)
 patchline(boundaries{ii}(:,2),boundaries{ii}(:,1),...
    'edgealpha',0.2,'edgecolor',cmap(ind(ii),:),'linewidth',3);
 end